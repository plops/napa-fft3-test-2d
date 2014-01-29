#.(require :ltk)
#.(require "napa-fft3")
(defpackage :l (:use :cl :ltk))
(in-package :l)


(defun read-pgm (filename)
  (declare ((or pathname string) filename)
           (values (or (simple-array (unsigned-byte 8) 2)
                       (simple-array (unsigned-byte 16) 2)) &optional))
  (with-open-file (s filename)
    (unless (equal (symbol-name (read s)) "P5")
      (error "no PGM file"))
    (read-line s)
    (let* ((w (read s))
           (h (read s))
           (grays (read s))
           (pos (file-position s)))
      (declare ((integer 0 65535) grays w h))
      (let* ((type (if (<= grays 255)
                       '(unsigned-byte 8)
                       '(unsigned-byte 16)))
             (data (make-array (list h w)
                               :element-type type))
             (data-1d (make-array (* h w)
                                  :element-type type
                                  :displaced-to data)))
        (with-open-file (s2 filename :element-type type)
          (file-position s2 pos)
          (read-sequence data-1d s2))
        data))))


(defun double (d)
  (let* ((n (reduce #'* (array-dimensions d)))
	 (a1 (make-array n
			 :element-type 'double-float))
	 (a (make-array (array-dimensions d)
			:element-type (array-element-type a1)
			:displaced-to a1))
	 (d1 (make-array n
			 :element-type (array-element-type d)
			 :displaced-to d)))
    (dotimes (i n)
      (setf (aref a1 i) (* (/ 255d0) (aref d1 i))))
    a))

(defun complex-double-float (d)
  (let* ((n (reduce #'* (array-dimensions d)))
	 (a1 (make-array n
			 :element-type '(complex double-float)))
	 (a (make-array (array-dimensions d)
			:element-type (array-element-type a1)
			:displaced-to a1))
	 (d1 (make-array n
			 :element-type (array-element-type d)
			 :displaced-to d)))
    (dotimes (i n)
      (setf (aref a1 i) (complex (* (/ 255d0) (aref d1 i)))))
    a))

(defun checker (a)
  "multiply with [-1 +1 -1 ...; +1 -1 +1; ....; ...] to shift fourier transform into middle"
  (let ((b (make-array (array-dimensions a)
		       :element-type (array-element-type a))))
   (destructuring-bind (h w) (array-dimensions a)
     (dotimes (i w)
       (dotimes (j h)
	 (setf (aref b j i) (* (expt -1 (+ i j)) (aref a j i))))))
   b))

(defun damp-edge (&key a (width .1))
  "force the values on the edge of the array to be zero"
  (destructuring-bind (h w) (array-dimensions a)
   (let ((b (make-array (array-dimensions a)
			:element-type (array-element-type a)))
	 (win-x (make-array w :element-type 'double-float
			    :initial-element 1d0))
	 (win-y (make-array h :element-type 'double-float
			    :initial-element 1d0)))
     (let ((l (floor (* width w))))
       (dotimes (i l)
	 (setf (aref win-x i) (sin (/ (* .5 pi i) l))))
       (dotimes (i l)
	 (setf (aref win-x (+ w -1 (- i))) (sin (/ (* .5 pi i) l)))))
     (let ((l (floor (* width h))))
       (dotimes (i l)
	 (setf (aref win-y i) (sin (/ (* .5 pi i) l))))
       (dotimes (i l)
	 (setf (aref win-y (+ h -1 (- i))) (sin (/ (* .5 pi i) l)))))
     (destructuring-bind (h w) (array-dimensions a)
       (dotimes (i w)
	 (dotimes (j h)
	   (setf (aref b j i) (* (aref win-x i) (aref win-y j) (aref a j i))))))
     b)))

(defun ubyte (d &key (scale 1d0))
  (let* ((n (reduce #'* (array-dimensions d)))
	 (a1 (make-array n
			 :element-type '(unsigned-byte 8)))
	 (a (make-array (array-dimensions d)
			:element-type (array-element-type a1)
			:displaced-to a1))
	 (d1 (make-array n
			 :element-type (array-element-type d)
			 :displaced-to d)))
    (dotimes (i n)
      (setf (aref a1 i) (min 255 (max 0 (floor (* scale (abs (aref d1 i))))))))
    a))



(defun write-pgm (filename img)
  (declare (simple-string filename)
           ((array (unsigned-byte 8) 2) img)
           (values null &optional))
  (destructuring-bind (h w)
      (array-dimensions img)
    (declare ((integer 0 65535) w h))
    (with-open-file (s filename
                       :direction :output
                       :if-exists :supersede
                       :if-does-not-exist :create)
      (declare (stream s))
      (format s "P5~%~D ~D~%255~%" w h))
    (with-open-file (s filename 
                       :element-type '(unsigned-byte 8)
                       :direction :output
                       :if-exists :append)
      (let ((data-1d (make-array 
                      (* h w)
                      :element-type '(unsigned-byte 8)
                      :displaced-to img)))
        (write-sequence data-1d s)))
    nil))

(defvar *dat* nil)

(defun fft2 ()
  "use 1d fft from napa-fft3 to transform a real 2d image. for now,
fft2 acts destructive on the data in the array *dat*"
  (destructuring-bind (h w) (array-dimensions *dat*)
   (let* ((a1 (make-array (* h w)
			 :element-type 'double-float
			 :displaced-to *dat*))
	 (outline (make-array w
			   :element-type '(complex double-float)))
	 (out1 (make-array (* w h)
			   :element-type '(complex double-float)))
	  (out2 (make-array (list w h)
			   :element-type '(complex double-float)
			   :displaced-to out1))
	  (outcol (make-array h
			   :element-type '(complex double-float)))
	  (out-final (make-array (list h w)
			    :element-type '(complex double-float)
			    )))
     (loop for j below h do
	  (let ((in (make-array w
				:element-type 'double-float
				:displaced-to a1
				:displaced-index-offset (* j w))))
	    (napa-fft:rfft in :size w :dst outline)
	    (dotimes (i w)
	      (setf (aref out2 i j) (aref outline i)))))
     (loop for i below w do
	  (let ((in (make-array h
				:element-type '(complex double-float)
				:displaced-to out1
				:displaced-index-offset (* i h))))
	    (napa-fft:fft in :size h :dst outcol)
	    (dotimes (j h)
	      (setf (aref out-final j i) (aref outcol j)))))
     (defparameter *out* out-final))))

(defun fft2c ()
  "use 1d fft from napa-fft3 to transform a complex 2d image. for now,
fft2 acts destructive on the data in the array *dat*"
  (destructuring-bind (h w) (array-dimensions *dat*)
   (let* ((a1 (make-array (* h w)
			 :element-type '(complex double-float)
			 :displaced-to *dat*))
	 (outline (make-array w
			   :element-type '(complex double-float)))
	  (out1 (make-array (* w h)
			    :element-type '(complex double-float)))
	  (out2 (make-array (list w h)
			    :element-type '(complex double-float)
			    :displaced-to out1))
	  (outcol (make-array h
			   :element-type '(complex double-float)))
	  (out-final (make-array (list h w)
			    :element-type '(complex double-float)
			    )))
     (loop for j below h do
	  (let ((in (make-array w
				:element-type '(complex double-float)
				:displaced-to a1
				:displaced-index-offset (* j w))))
	    (napa-fft:fft in :size w :dst outline)
	    (dotimes (i w)
	      (setf (aref out2 i j) (aref outline i)))))
     (loop for i below w do
	  (let ((in (make-array h
				:element-type '(complex double-float)
				:displaced-to out1
				:displaced-index-offset (* i h))))
	    (napa-fft:fft in :size h :dst outcol)
	    (dotimes (j h)
	      (setf (aref out-final j i) (aref outcol j)))))
     (defparameter *out* out-final))))

(declaim (optimize (speed 0) (safety 3) (debug 3)))

(defun gauss-blur2c (&key a (sigma 15d0))
  "gaussian convolution of a 2d array. convolve each line and then
each column by multiplying with a window in the fourier domain. no
time-consuming bit-reversals are necessary.
A .. 2d input image, its contents will be destroyed and replaced with
     the result.
SIGMA .. width of the gaussian in pixels"
  (declare (type (array (complex double-float) 2) a))
  (destructuring-bind (h w) (array-dimensions a)
    (let ((winx (make-array w :element-type '(complex double-float)))
	  (row (make-array w :element-type '(complex double-float)))
	  (row-out (make-array w :element-type '(complex double-float))))
      (let ((n (/ (loop for i below w sum (exp (/ (- (expt i 2))
						  (expt sigma 2)))))))
       (dotimes (i w)
	 (setf (aref winx i) (complex (* n (exp (/ (- (expt i 2))
						    (expt sigma 2))))))))
     
     (let ((kwinx (napa-fft:fft winx :in-order nil)))
       (dotimes (j h)
	 (let ((a1 (make-array w :element-type '(complex double-float)
			       :displaced-to a :displaced-index-offset (* j w))))
	   (napa-fft:fft a1 :in-order nil :dst row)
	   (napa-fft:ifft row :in-order nil :window kwinx :dst row-out)
	   (dotimes (i w) ;; unfortunately, ifft wants a simple-array as dst
	     (setf (aref a j i) (aref row-out i)))))))
    
    (let ((winy (make-array h :element-type '(complex double-float)))
	  (col (make-array h :element-type '(complex double-float)))
	  (col-in (make-array h :element-type '(complex double-float))))
      (let ((n (/ (loop for j below h sum (exp (/ (- (expt j 2))
						  (expt sigma 2)))))))
	(dotimes (j h)
	  (setf (aref winy j) (complex (* n (exp (/ (- (expt j 2))
						    (expt sigma 2))))))))
      (let ((kwiny (napa-fft:fft winy :in-order nil)))
	(dotimes (i w)
	  (dotimes (j h)
	    (setf (aref col-in j) (aref a j i)))
	  (napa-fft:fft col-in :in-order nil :dst col)
	  (napa-fft:ifft col :in-order nil :window kwiny :dst col-in)
	  (dotimes (j h)
	    (setf (aref a j i) (aref col-in j)))))
      a)))

(defun .linear (a)
  (make-array (reduce #'* (array-dimensions a))
	      :element-type (array-element-type a)
	      :displaced-to a))

(defun .* (a b)
  (let* ((c (make-array (array-dimensions a) :element-type (array-element-type a)))
	(a1 (.linear a))
	(b1 (.linear b))
	(c1 (.linear c)))
    (dotimes (i (length c))
      (setf (aref c1 i) (* (aref a1 i) (aref b1 i))))
    c))

#+nil
(time
 (let ((files (directory "/home/martin/dat/r/*.pgm")))
   (let ((e (elt files 327)))
					;loop for e in files do
       
     (let ((base (string-trim (list #\/) (pathname-name e))))
       (defparameter *dat* (extract :a (complex-double-float (read-pgm e))))
       (write-pgm (format nil "/dev/shm/f.pgm") (ubyte *dat* :scale 255d0))
       (write-pgm (format nil "/dev/shm/f2.pgm") (ubyte (gauss-blur2c :a *dat*) :scale 255d0))
       ))))


(defun next-power-of-two (n)
  (expt 2 (ceiling (log n 2))))

(defun extract (&key
		  a
                (x (floor (array-dimension a 1) 2))
                (y (floor (array-dimension a 0) 2)) 
                (w (next-power-of-two (array-dimension a 1) ;(max x (- (array-dimension a 1) x))
		    ))
                (h (next-power-of-two (array-dimension a 0) ;(max y (- (array-dimension a 0) y))
				      )))
  (destructuring-bind (hh ww) (array-dimensions a)
   (let* ((zero (coerce 0 (array-element-type a)))
	  (b1 (make-array (* h w) :element-type (array-element-type a)
			  :initial-element zero))
	  (b (make-array (list h w)
			 :element-type (array-element-type a)
			 :displaced-to b1))
	  (ox (- x (floor w 2)))
	  (oy (- y (floor h 2))))
     ;; (assert (<= 0 ox))
     ;; (assert (<= 0 oy))
     ;; (assert (< (+ w ox) (array-dimension a 1)))
     ;; (assert (< (+ h oy) (array-dimension a 0)))
     (dotimes (j h)
       (dotimes (i w)
	 (setf (aref b j i)
	       (let ((ii (+ i ox))
		     (jj (+ j oy)))
		 (if (and (<= 0 jj (1- hh)) (<= 0 ii (1- ww)))
		     (aref a jj ii)
		     zero)))))
     b)))


#+nil
(progn (defparameter *dat* (damp-edge :a (checker (double (read-pgm (first
								     (directory "~/dat/bla0*.pgm")))))))
       (write-pgm "/dev/shm/o2.pgm" (ubyte *dat* :scale 255d0))
       (fft2)
       
       
       (write-pgm "/dev/shm/o.pgm" (ubyte *out* :scale 1d0))) 

#+nil
(let ((files (directory "/home/martin/dat/r/*.pgm")))
  (; let ((e (elt files 327)))
   loop for e in files do
       
       (let ((base (string-trim (list #\/) (pathname-name e))))
	 (defparameter *dat* (extract :a (checker (damp-edge :width .1 :a (double (read-pgm e))))))
	 (write-pgm (format nil "/dev/shm/f.pgm") (ubyte *dat* :scale 255d0))
	 (fft2)
	 (write-pgm (format nil "/dev/shm/o.pgm") (ubyte *out* :scale 1d0))
	 (let ((small (damp-edge :width .1 :a (extract :x 370 :y 490 :w 128 :h 128 :a *out*))))
	   (write-pgm (format nil "/dev/shm/o~a.pgm" base) (ubyte small :scale 1d0))
	   (setf *dat* small) (fft2c)
	   (write-pgm (format nil "/dev/shm/y~a.pgm" base) (ubyte *out* :scale .01d0))))))
;; make mosaic from the reconstructed images
;; montage  y* -tile 20x18 -geometry 128x128+0+0  mon.png


;; roi 121x89+546+390
#+nil
(list (+ 546 (floor 121 2))
      (+ 390 (floor 89 2)))
#+nil
(list (/ (- (+ 546 (floor 121 2)) 256)
	 1d0)
      (/ (- (+ 390 (floor 89 2)) 512)
	 1d0))

#+nil
(progn (defparameter *dat* (extract (checker (double (read-pgm (elt
								
								59))))))
       (write-pgm "/dev/shm/o2.pgm" (ubyte *dat* :scale 255d0))
       (fft2)
       
       
       (write-pgm "/dev/shm/o.pgm" (ubyte *out* :scale 1d0)))

(defun read-and-discard-comments (s)
  (let ((l (read-line s)))
    (unless (eq #\# (aref l 0))
      l)))

(read-char)

(defparameter *eof-object* (make-symbol "EOF-OBJECT"))
(defun eofp (char)
  (eq char *eof-object*))
(defun token-delimiterp (char)
  (member char (list #\Newline #\Space #\Tab #\Return)))
(defun token-commentp (char)
  (eq char #\#))

(defun digitp (char)
  (member char (list #\0 #\1 #\2 #\3 #\4 #\5 #\6 #\7 #\8 #\9)))

(defun read-separate-and-process-whitespace (s)
  (let ((res nil))
    (loop for c = (read-char s nil *eof-object*) while
	 (cond ((eofp c) nil)
	       ((token-delimiterp c) (unread-char c s) t)
	       ((token-commentp c) (unread-char c s) nil)
	       (t t))
	 do
	 (push c res))
    (reverse res)))

(defun read-token (s)
  (let ((c (read-char s nil *eof-object*)))
    (unread-char c s)
    (list (cond ((digitp c) :number)
		((eq #\P c) :magic)
		(t :unknown))
	  (read-separate-and-process-whitespace s))))

(with-open-file (s  (elt
		     (directory "/dev/shm/crop/*.pgm")
		     32))
  (list   (read-token s)
	  (read-token s)
	  (read-token s)
	  (read-token s)))

(defvar *out* nil)
(defun canvas-test ()
  (with-ltk ()
    (let* ((c (make-instance 'canvas :width (+ 515 512 ) :height 512
			     :cursor "crosshair"
			     ))
	   
	   (dat2 (image-load (make-image)
			    (first
			     (directory "/dev/shm/o2.pgm"))))
	   (dat (image-load (make-image)
			    (first
			     (directory "/dev/shm/o.pgm"))))
	   (im (create-image c 0 0 :image dat))
	   (r (let* ((x 286)
		    (y 54)
		    (w 76)
		    (h w))
		(create-rectangle c (- x (floor w 2))
				  (- y (floor h 2))
				  (+ x (floor w 2))
				  (+ y (floor h 2)))))
	   (im2 (create-image c 515 0 :image dat2))
	   (text (create-text c 10 10 "10 10")))
      
      (bind c "<Motion>"
	    #'(lambda (evt)
		(itemdelete c text)
		(setf text
		      (create-text
		       c
		       (+ 4 (event-x evt))
		       (+ 3 (event-y evt))
		       (format nil "~4d ~4d ~a~%"
			       (event-x evt) (event-y evt)
			       (when *out*
				 (destructuring-bind (h w)
				     (array-dimensions *out*)
				   (let ((z 
					  (aref *out*
						(min (1- h)
						     (max 0
							  (event-y evt)))
						(min (1- w)
						     (max 0
							  (event-x evt))))))
				     (format nil "~6,1f ~6,1f"
					     (abs z)
					     (* 180 (/ pi) (phase z)))))))))
			     (itemconfigure c text 'fill "red")))
      (itemconfigure c r 'outline "red")
      (pack c :expand 1 :fill :both))))

#+nil
(canvas-Test)

#+nil
(ltk::ltktest)

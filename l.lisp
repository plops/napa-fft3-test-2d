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

(defun checker (a)
  "multiply with [-1 +1 -1 ...; +1 -1 +1; ....; ...] to shift fourier transform into middle"
  (let ((b (make-array (array-dimensions a)
		       :element-type (array-element-type a))))
   (destructuring-bind (h w) (array-dimensions a)
     (dotimes (i w)
       (dotimes (j h)
	 (setf (aref b j i) (* (expt -1 (+ i j)) (aref a j i))))))
   b))

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



(defun fft2 ()
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

#+nil
(progn (defparameter *dat* (checker (double (read-pgm (first
					       (directory "~/dat/bla0*.pgm"))))))
       (fft2)
       

       (write-pgm "/dev/shm/o.pgm" (ubyte *out* :scale 1d0))) 


(defun canvas-test ()
  (with-ltk ()
    (let* ((c (make-instance 'canvas))
	   (line (create-line c (list 100 100 400 50 700 165)))
	   (dat (image-load (make-image)
			    (first
			     (directory "~/dat/bla0*.pgm"))))
	   (dat2 (image-load (make-image)
			    (first
			     (directory "/dev/shm/o.pgm"))))
	   (im (create-image c 0 0 :image dat))
	   (im2 (create-image c 512 0 :image dat2)))
      (pack c :expand 1 :fill :both))))

#+nil
(canvas-Test)

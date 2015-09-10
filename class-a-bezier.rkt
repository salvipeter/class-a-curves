;;; Class A Bezier Curve Test
;;; by Peter Salvi, 2015.09.10.

;;; Based on the paper by Yoshida et al. (2008)

#lang racket

;;; Libraries
(require racket/gui)

;;; Parameters
(define point-radius 4)
(define line-width 0)
(define resolution 100)

;;; Default placements
(define a0 '(44 340))
(define a1 '(235 88))
(define a2 '(586 279))

;;; Options
(define optimization-iterations 50)
(define epsilon 1e-8)

;;; Variables
(define degree 4)
(define dragged #f)
(define show-bezier-cpts? #t)
(define bezier-cpts '())

;;; Basic Maths
(define (binomial n k)
  (if (= k 0)
      1
      (* (/ n k) (binomial (- n 1) (- k 1)))))
(define (v+ . args) (apply map + args))
(define (v- . args) (apply map - args))
(define (v* u . args) (map (lambda (x) (apply * x args)) u))
(define (vlength u) (sqrt (apply + (map (lambda (x) (* x x)) u))))
(define (vnormalize u) (v* u (/ (vlength u))))
(define (point-distance p q) (vlength (v- q p)))
(define (scalar-product u v) (apply + (map * u v)))
(define (cross-product u v)
  (list (- (* (second u) (third v)) (* (third u) (second v)))
	(- (* (third u) (first v)) (* (first u) (third v)))
	(- (* (first u) (second v)) (* (second u) (first v)))))
(define (to-system u d)
  (let* ([u (vnormalize u)]
         [v (list (- (second u)) (first u))])
    (list (scalar-product d u) (scalar-product d v))))
(define (from-system u d)
  (let* ([u (vnormalize u)]
         [v (list (- (second u)) (first u))])
    (v+ (v* u (first d)) (v* v (second d)))))

;;; Curves

(define (bezier-eval-one-point points u)
  (let* ([n (- (length points) 1)]
         [v (- 1 u)]
         [p '(0 0)])
    (for ([k (in-range (+ n 1))]
          [q points])
      (set! p (v+ p (v* q (binomial n k) (expt u k) (expt v (- n k))))))
    p))

(define (bezier-evaluate points)
  (for/list ([i (in-range resolution)])
    (let* ([u (/ i (- resolution 1))]
           [p (bezier-eval-one-point points u)])
      (make-object point% (first p) (second p)))))

(define (rotate p alpha)
  (let ([c (cos alpha)]
        [s (sin alpha)]
        [x (first p)]
        [y (second p)])
    (list (- (* c x) (* s y))
          (+ (* s x) (* c y)))))

(define (angle-between u v)
  (* (acos (min 1 (max -1 (scalar-product (vnormalize u)
                                          (vnormalize v)))))
     (let ([cross (cross-product (append u '(0)) (append v '(0)))])
       (if (< (scalar-product '(0 0 1) cross) 0)
         -1
         1))))

(define (evaluate-polynomial as x)
  (if (empty? as)
      0
      (+ (first as)
         (* x (evaluate-polynomial (rest as) x)))))

(define (multiply-all-but as x)
  (if (empty? as)
      1
      (* (if (= (first as) x)
             1 
             (- x (first as)))
         (multiply-all-but (rest as) x))))

(define (find-polynomial-roots as)
  "Given a list of numbers (a0 a1 a2 ... an-1),
find x such that sum ai * x^i = 0.
Uses the Durand-Kerner method."
  (let* ([n (- (length as) 1)]
         [xs (for/list ([i (range n)])
               (expt 0.4+0.9i i))]) ; random complex number
    (for ([i (range optimization-iterations)])
      (set! xs (map (lambda (x)
                      (- x (/ (evaluate-polynomial as x)
                              (multiply-all-but xs x))))
                    xs)))
    xs))

;;; for example: (find-polynomial-roots '(-5 3 -3 1))

(define (closest-to-one lst)
  (first (sort lst (lambda (x y)
                     (< (abs (- (abs x) 1))
                        (abs (- (abs y) 1)))))))

(define (find-polynomial-root as)
  "Choose one root, that is real and positive."
  (let* ([xs (find-polynomial-roots as)]
         [ys (filter (lambda (x) (< (imag-part x) epsilon)) xs)]
         [zs (map real-part ys)]
         [ws (filter (lambda (z) (> z 0)) zs)])
    (if (empty? ws)
        (if (empty? zs)
            1 ; kutykurutty
            (closest-to-one zs))
        (closest-to-one ws))))

(define (other-end coeffs s)
  (define (rec vs i)
    (if (empty? vs)
        '(0 0)
        (v+ (v* (first vs) (expt s i))
            (rec (rest vs) (+ i 1)))))
  (rec coeffs 0))

(define (optimize coeffs target)
  "Given COEFFS as a list of (v0 v1 v2 ... vn-1),
find s such that sum s^i * vi has the same direction as TARGET.
The return value is (s b0), such that b0 * sum s^i * vi = TARGET.
Note that all coefficients and the target are two-dimensional."
  (let* ([vs (map (lambda (p)
                    (second (to-system target p)))
                  coeffs)]
         [s (find-polynomial-root vs)])
    (list s (/ (vlength target)
               (vlength (other-end coeffs s))))))

(define (generate-bezier-cpts v s alpha)
  "Generates the Bezier control points.
Uses also A0 and DEGREE."
  (define (rec i p u)
    (if (> i degree)
        '()
        (cons p (rec (+ i 1) (v+ p u) (v* (rotate u alpha) s)))))
  (rec 0 a0 v))

(define (generate-rotations u alpha k)
  (if (= k 1)
      (list u)
      (cons u (generate-rotations (rotate u alpha) alpha (- k 1)))))

(define (recompute-curve)
  "Generate BEZIER-CPTS based on A0, A1, A2 and DEGREE."
  (let* ([alpha (/ (angle-between (v- a1 a0) (v- a2 a1)) (- degree 1))]
         [u (vnormalize (v- a1 a0))]
         [vs (generate-rotations u alpha degree)]
         [s-b0 (optimize vs (v- a2 a0))]
         [s (first s-b0)]
         [b0 (second s-b0)]
         [v (v* u b0)])
    (set! bezier-cpts (generate-bezier-cpts v s alpha))))

;;; Graphics

(define (draw-point dc p)
  (send dc draw-ellipse
        (- (first p) point-radius) (- (second p) point-radius)
        (* point-radius 2) (* point-radius 2)))

(define (draw-segment dc p q)
  (send dc draw-line (first p) (second p) (first q) (second q)))

(define (draw canvas dc)
  (send dc set-pen "GREEN" line-width 'solid)
  (send dc draw-lines (bezier-evaluate bezier-cpts))
  (send dc set-pen "BLUE" line-width 'solid)
  (draw-segment dc a0 a1) (draw-segment dc a1 a2)
  (unless (or (empty? bezier-cpts) (not show-bezier-cpts?))
    (send dc set-brush "RED" 'solid)
    (send dc set-pen "RED" line-width 'solid)
    (for ([p bezier-cpts] [q (rest bezier-cpts)])
      (draw-segment dc p q))
    (for-each (lambda (p) (draw-point dc p)) bezier-cpts))
  (send dc set-brush "BLACK" 'solid)
  (send dc set-pen "BLACK" line-width 'solid)
  (for-each (lambda (p) (draw-point dc p)) (list a0 a1 a2)))

;;; GUI

(define (handle-mouse-movement event)
  (if dragged
      (let ([p (list (send event get-x) (send event get-y))])
        (case dragged
          [(0) (set! a0 p) (recompute-curve)]
          [(1) (set! a1 p) (recompute-curve)]
          [(2) (set! a2 p) (recompute-curve)])
        #t)
    #f))

(define (handle-mouse-down event)
  (if dragged
      (handle-mouse-up event)
      (let ([p (list (send event get-x) (send event get-y))]
            [points (list a0 a1 a2)])
        (for ([i (in-range 3)]
              [pi points])
          (when (< (point-distance p pi) point-radius)
            (set! dragged i))))))

(define (handle-mouse-up event)
  (set! dragged #f)
  #t)

(define scurve-canvas%
  (class canvas%
    (inherit refresh)
    (define/override (on-event event)
      (when (case (send event get-event-type)
              [(motion) (handle-mouse-movement event)]
              [(left-down) (handle-mouse-down event)]
              [(left-up) (handle-mouse-up event)])
        (refresh)))
    (super-new)))

(let* ([frame (new frame% [label "Class A Bezier Curves"])]
       [vbox (new vertical-pane% [parent frame])]
       [canvas (new scurve-canvas% [parent vbox]
                    [min-width 640] [min-height 480]
                    [paint-callback draw])]
       [hbox1 (new horizontal-pane% [parent vbox])]
       [label (new message% [label ""] [parent hbox1]
                   [auto-resize #t])]
       [set-label (lambda ()
                    (send label set-label
                          (string-append "Degree: "
                                         (number->string degree))))]
       [slider-tmp #f]
       [slider (new slider% [label ""] [min-value 3] [max-value 12] [parent hbox1]
                    [init-value degree] [style '(horizontal plain)]
                    [min-width 200] [stretchable-width #f]
                    [callback (lambda (c e)
                                (set! degree (send slider-tmp get-value))
                                (set-label)
                                (recompute-curve)
                                (send canvas refresh))])])
  (let-syntax ([add-option (syntax-rules ()
			     [(_ var a-label box)
			      (new check-box% [parent box]
                                   [label a-label] [value var]
				   [callback (lambda (c e)
					       (set! var (not var))
					       (send canvas refresh))])])])
    (add-option show-bezier-cpts? "Bezier control points" hbox1))
  (set! slider-tmp slider)
  (set-label)
  (recompute-curve)
  (send frame show #t))

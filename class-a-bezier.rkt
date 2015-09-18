;;; Class A Curve Test
;;; by Peter Salvi, September 2015.

;;; Bezier curves based on the paper by Yoshida et al. (2008)
;;; LA curves based on the paper by Miura (2006) and Yoshida (2006)

#lang racket

;;; Libraries
(require math/matrix)
(require racket/gui)

;;; Parameters
(define point-radius 4)
(define line-width 0)
(define resolution 100)
(define optimization-iterations 100)
(define bisection-iterations 10)
(define epsilon 1e-5)

;;; Default placements
(define a0 '(44 340))
(define a1 '(235 88))
(define a2 '(586 279))

;;; Variables
(define dragged #f)
(define degree 4)
(define show-bezier-cpts? #t)
(define alpha -1)
(define show-lac? #f)

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
(define (to-system u d)
  (let* ([u (vnormalize u)]
         [v (list (- (second u)) (first u))])
    (list (scalar-product d u) (scalar-product d v))))
(define (from-system u d)
  (let* ([u (vnormalize u)]
         [v (list (- (second u)) (first u))])
    (v+ (v* u (first d)) (v* v (second d)))))

;;; Bezier Curve

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

(define (clockwise? u v)
  (< (- (* (first u) (second v)) (* (second u) (first v)))
     0))

(define (angle-between u v)
  (* (acos (min 1 (max -1 (scalar-product (vnormalize u)
                                          (vnormalize v)))))
     (if (clockwise? u v) -1 1)))

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
  (if (< (last as) 0)
      (find-polynomial-roots (map - as))
      (let* ([n (- (length as) 1)]
             [xs (for/list ([i (range n)])
                   (expt 0.4+0.9i i))]) ; random complex number
        (for ([i (range optimization-iterations)])
          (set! xs (map (lambda (x)
                          (- x (/ (evaluate-polynomial as x)
                                  (multiply-all-but xs x))))
                        xs)))
        xs)))

(define (minimal-imag-part lst)
  (real-part
   (first (sort lst
                (lambda (x y)
                  (< (abs (imag-part x))
                     (abs (imag-part y))))))))

(define (closest-to-zero as lst)
  (first (sort lst (lambda (x y)
                     (< (abs (evaluate-polynomial as x))
                        (abs (evaluate-polynomial as x)))))))

(define (find-polynomial-root as)
  "Choose one root, that is real and positive."
  (let* ([xs (find-polynomial-roots as)]
         [ys (filter (lambda (x) (< (imag-part x) epsilon)) xs)]
         [zs (map real-part ys)]
         [ws (filter (lambda (z) (> z 0)) zs)])
    (if (empty? ws)
        (minimal-imag-part xs) ; kutykurutty
        (closest-to-zero as ws))))

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
    (if (= i degree)
        (list p)
        (cons p (rec (+ i 1) (v+ p u) (v* (rotate u alpha) s)))))
  (rec 0 a0 v))

(define (generate-rotations u alpha k)
  (if (= k 1)
      (list u)
      (cons u (generate-rotations (rotate u alpha) alpha (- k 1)))))

(define (class-a-bezier-fit)
  "Generate control points based on A0, A1, A2 and DEGREE."
  (let* ([angle (angle-between (v- a1 a0) (v- a2 a1))]
         [alpha (/ angle (- degree 1))]
         [u (vnormalize (v- a1 a0))]
         [vs (generate-rotations u alpha degree)]
         [s-b0 (optimize vs (v- a2 a0))]
         [s (first s-b0)]
         [b0 (second s-b0)]
         [v (v* u b0)])
    (generate-bezier-cpts v s alpha)))

;;; LA curve

(define gaussian-quadrature
  '((-0.861136312 0.347854845) (-0.339981044 0.652145155)
    (0.339981044 0.652145155) (0.861136312 0.347854845)))

(define (integrate fn a b)
  (let ([c1 (/ (- b a) 2)]
        [c2 (/ (+ a b) 2)])
    (* c1 (apply + (map (lambda (x-w)
                          (* (second x-w)
                             (fn (+ (* (first x-w) c1) c2))))
                        gaussian-quadrature)))))

(define (complex->point p)
  (list (real-part p) (imag-part p)))

(define (complex->canvas-point p)
  (make-object point% (real-part p) (imag-part p)))

(define (point->complex p)
  (+ (first p) (* 0+1i (second p))))

(define (compute-theta-ef l theta-d)
  "Compute THETA-E / THETA-F given Lambda of a Standard Form II curve."
  (let ([tri (compute-triangle l theta-d)])
    (angle-between (complex->point (- (second tri) (first tri)))
                   (complex->point (- (third tri) (first tri))))))

(define (bisection theta-d theta-e theta-f)
  "As in Appendix A of Yoshida'06."
  (let ([enlarge (= alpha 1)])
    (define (rec lmin lmax i)
      (let ([l (/ (+ lmin lmax) 2)])
        (if (= i 0)
            l ; not found - return best estimate
            (let ([f (if (<= alpha 1)
                         (- theta-e (compute-theta-ef l theta-d))
                         (- theta-f (compute-theta-ef l theta-d)))])
              (cond [(< (abs f) epsilon) l]
                    [(or (and (> f 0) (<= alpha 1))
                         (and (< f 0) (> alpha 1)))
                     (rec l (if enlarge (* lmax 10) lmax) (- i 1))]
                    [else
                     (set! enlarge #f)
                     (rec lmin l (- i 1))])))))
    (let ([lmax (cond [(< alpha 1) (/ (* theta-d (- 1 alpha)))]
                      [(> alpha 1) (/ (* -1 theta-d (- 1 alpha)))]
                      [else 1])])
      (rec 0 lmax bisection-iterations))))

(define (map-points tri)
  "Returns a mapping that maps this triangle to A0,A1,A2."
  (let* ([b0 (list (real-part (first tri)) (imag-part (first tri)))]
         [b1 (list (real-part (second tri)) (imag-part (second tri)))]
         [b2 (list (real-part (third tri)) (imag-part (third tri)))]
         [M (list*->matrix (list (list (first b0) (first b1) (first b2))
                                 (list (second b0) (second b1) (second b2))
                                 '(1 1 1)))]
         [A (matrix* (list*->matrix (list (list (first a0) (first a1) (first a2))
                                          (list (second a0) (second a1) (second a2))
                                          '(1 1 1)))
                     (matrix-inverse M))])
    (lambda (p)
      (let* ([v (->col-matrix (list (real-part p) (imag-part p) 1))]
             [q (matrix* A v)])
        (+ (matrix-ref q 0 0) (* 0+1i (matrix-ref q 1 0)))))))

(define (theta-d->s l theta)
  (case alpha
    [(0) (* (- l) (log (- 1 theta)))]
    [(1) (/ (- (exp (* theta l)) 1) l)]
    [else (/ (- (expt (+ 1 (* theta l (- alpha 1)))
                      (/ alpha (- alpha 1)))
                1)
             (* l alpha))]))

(define (compute-triangle l theta-d)
  (let* ([p0 0]
         [p2 (lac-eval-one-point-yoshida l 0 (theta-d->s l theta-d))]
         [eith (exp (* 0+1i theta-d))]
         [gamma (- (/ (imag-part p2) (imag-part eith)))]
         [p1 (real-part (+ p2 (* eith gamma)))]) ; intersection of P0P1 and P1P2
    (list p0 p1 p2)))

(define (lac-fit)
  "Returns (L MAPPING S-FROM S-TO)."
  (let* ([theta-d (angle-between (v- a1 a0) (v- a2 a1))]
         [theta-e (angle-between (v- a1 a0) (v- a2 a0))]
         [theta-f (angle-between (v- a1 a2) (v- a0 a2))]
         [l (bisection theta-d theta-e theta-f)]
         [tri (compute-triangle l theta-d)])
    (if (> alpha 1)
        (list l (map-points (reverse tri)) (theta-d->s l (- theta-d)) 0)
        (list l (map-points tri) 0 (theta-d->s l theta-d)))))

(define (lac-eval-one-point-miura params s-from s-to)
  "Does not work in this form for alpha=0,1. Not used now."
  (let* ([p0 (first params)]
         [c0 (second params)]
         [c1 (third params)]
         [c2 (fourth params)]
         [theta (lambda (s)
                  (+ (/ (* alpha (expt (+ (* c0 s) c1)
                                       (/ (- alpha 1) alpha)))
                        (* (- alpha 1) c0))
                     c2))]
         [integrand (lambda (s) (exp (* 0+1i (theta s))))])
    (+ p0 (integrate integrand s-from s-to))))

(define (lac-eval-one-point-yoshida l s-from s-to)
  (let* ([theta (lambda (s)
                  (case alpha
                    [(0) (- 1 (exp (* -1 l s)))]
                    [(1) (/ (log (+ (* l s) 1)) l)]
                    [else (/ (- (expt (+ (* l alpha s) 1)
                                      (- 1 (/ alpha)))
                                1)
                             (* l (- alpha 1)))]))]
         [integrand (lambda (s) (exp (* 0+1i (theta s))))])
    (integrate integrand s-from s-to)))

(define (lac-evaluate params)
  "PARAMS is a list of (L MAPPING S-FROM S-TO)."
  (let ([pa (point->complex a0)]
        [pb (point->complex a1)]
        [pc (point->complex a2)])
    (when (> (magnitude (* pa pb)) (magnitude (* pb pc)))
      (set! a0 (complex->point pc))
      (set! a2 (complex->point pa))))
  (let* ([l (first params)]
         [mapping (second params)]
         [s-from (third params)]
         [s-to (fourth params)]
         [s (for/list ([i (range resolution)])
              (+ s-from (* (- s-to s-from) (/ i (- resolution 1)))))]
         [points (map (lambda (x) (lac-eval-one-point-yoshida l s-from x)) s)])
    (map complex->canvas-point (map mapping points))))

;;; Graphics

(define (draw-point dc p)
  (send dc draw-ellipse
        (- (first p) point-radius) (- (second p) point-radius)
        (* point-radius 2) (* point-radius 2)))

(define (draw-segment dc p q)
  (send dc draw-line (first p) (second p) (first q) (second q)))

(define (draw canvas dc)
  (let ([bezier-cpts (class-a-bezier-fit)])
    (send dc set-pen "GREEN" line-width 'solid)
    (send dc draw-lines (bezier-evaluate bezier-cpts))
    (send dc set-pen "BLUE" line-width 'solid)
    (draw-segment dc a0 a1) (draw-segment dc a1 a2)
    (when show-bezier-cpts?
      (send dc set-brush "RED" 'solid)
      (send dc set-pen "RED" line-width 'solid)
      (for ([p bezier-cpts] [q (rest bezier-cpts)])
        (draw-segment dc p q))
      (for-each (lambda (p) (draw-point dc p)) bezier-cpts)))
  (when show-lac?
    (send dc set-pen "MAGENTA" line-width 'solid)
    (send dc draw-lines (lac-evaluate (lac-fit))))
  (send dc set-brush "BLACK" 'solid)
  (send dc set-pen "BLACK" line-width 'solid)
  (for-each (lambda (p) (draw-point dc p)) (list a0 a1 a2)))

;;; GUI

(define (handle-mouse-movement event)
  (if dragged
      (let ([p (list (send event get-x) (send event get-y))])
        (case dragged
          [(0) (set! a0 p)]
          [(1) (set! a1 p)]
          [(2) (set! a2 p)])
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
       [hbox2 (new horizontal-pane% [parent vbox])])
  (let-syntax ([add-option (syntax-rules ()
			     [(_ var a-label box)
			      (new check-box% [parent box]
                                   [label a-label] [value var]
				   [callback (lambda (c e)
					       (set! var (not var))
					       (send canvas refresh))])])]
               [add-slider (syntax-rules ()
                             [(_ var a-label box a-min a-max from-int to-int)
                              (let* ([label (new message% [label ""] [parent box]
                                                 [auto-resize #t])]
                                     [set-label (lambda ()
                                                  (send label set-label
                                                        (string-append a-label
                                                                       (number->string var))))]
                                     [slider-tmp #f]
                                     [slider (new slider% [label ""] [parent box]
                                                  [min-value (to-int a-min)]
                                                  [max-value (to-int a-max)]
                                                  [init-value (to-int var)]
                                                  [style '(horizontal plain)]
                                                  [min-width 200] [stretchable-width #f]
                                                  [callback (lambda (c e)
                                                              (let ([v (send slider-tmp get-value)])
                                                                (set! var (+ 0.0 (from-int v))))
                                                              (set-label)
                                                              (send canvas refresh))])])
                                (set! slider-tmp slider)
                                (set-label))])])
    (add-slider degree "Degree: " hbox1 3 12 identity identity)
    (add-option show-bezier-cpts? "Bezier control points" hbox1)
    (add-slider alpha "Alpha: " hbox2 -1 2 (lambda (x) (/ x 10)) (lambda (x) (* x 10)))
    (add-option show-lac? "Log-aesthetic curve" hbox2))
  (send frame show #t))

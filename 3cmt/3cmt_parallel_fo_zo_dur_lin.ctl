;; ==========================================================================
;; 3CMT_PARALLEL_FO_ZO_DUR_LIN
;; #34 of the three-compartment set
;;   Parallel absorption with first-order process and zero-order with duration
;;   Linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   DF x Dose     --(zero-order, D2)--> [Central (V2)]
;;   (1-DF) x Dose --> [Depot] --Ka--^
;;   [Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
;;   elimination from central: CL
;;
;; Data set / dosing requirements
;;   TWO dose records per administration, both with the full AMT:
;;     CMT=2 (central), RATE=-2  -> DF     x Dose, zero-order over D2
;;     CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, absorbed with Ka
;;                                  (starts immediately - parallel)
;;
;; Solver: closed-form (analytical)
;; ==========================================================================

$PROBLEM 3CMT parallel_fo_zo_dur | Linear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN12 TRANS4

; compartment 1 = depot
; compartment 2 = central (observation)
; compartment 3 = peripheral 1
; compartment 4 = peripheral 2
; compartment 5 = output

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V2    = THETA(2)*EXP(ETA(2))
  Q3    = THETA(3)*EXP(ETA(3))
  V3    = THETA(4)*EXP(ETA(4))
  Q4    = THETA(5)*EXP(ETA(5))
  V4    = THETA(6)*EXP(ETA(6))
  KA    = THETA(7)*EXP(ETA(7))
  D2    = THETA(8)*EXP(ETA(8))
  LGTDF = LOG(THETA(9)/(1-THETA(9))) + ETA(9)
  DF    = EXP(LGTDF)/(1+EXP(LGTDF))

; ---- structural / input specification ---------------------------------
  F1 = 1 - DF
  F2 = DF

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V2

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- initial estimates -------------------------------------------------
$THETA  (0, 5.0)             ; 1 CL     clearance (L/h)
$THETA  (0, 50.0)            ; 2 V2     central volume (L)
$THETA  (0, 10.0)            ; 3 Q3     inter-cmt clearance, periph.1 (L/h)
$THETA  (0, 100.0)           ; 4 V3     peripheral volume 1 (L)
$THETA  (0, 2.0)             ; 5 Q4     inter-cmt clearance, periph.2 (L/h)
$THETA  (0, 200.0)           ; 6 V4     peripheral volume 2 (L)
$THETA  (0, 1.0)             ; 7 KA     absorption rate constant (1/h)
$THETA  (0, 2.0)             ; 8 D2     zero-order duration into cmt 2 (h)
$THETA  (0.001, 0.5, 0.999)  ; 9 DF     dose fraction (logit-scale IIV)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V2
$OMEGA  0.09             ; 3 IIV Q3
$OMEGA  0.09             ; 4 IIV V3
$OMEGA  0.09             ; 5 IIV Q4
$OMEGA  0.09             ; 6 IIV V4
$OMEGA  0.09             ; 7 IIV KA
$OMEGA  0 FIX            ; 8 IIV D2    (free it if the data support it)
$OMEGA  0 FIX            ; 9 IIV DF    (free it if the data support it)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=3cmt_parallel_fo_zo_dur_lin.tab
$TABLE   ID CL V2 Q3 V3 Q4 V4 KA D2 DF ETA1 ETA2 ETA3 ETA4 ETA5 ETA6
         ETA7 ETA8 ETA9
         FIRSTONLY ONEHEADER NOPRINT FILE=3cmt_parallel_fo_zo_dur_lin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=3cmt_parallel_fo_zo_dur_lin.tab


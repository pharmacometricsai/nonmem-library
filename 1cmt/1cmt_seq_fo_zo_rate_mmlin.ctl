;; ==========================================================================
;; 1CMT_SEQ_FO_ZO_RATE_MMLIN
;; #33 of the one-compartment set
;;   Sequential absorption with first-order process with lag time followed by
;;     zero-order with rate
;;   Nonlinear and linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   DF x Dose     --(zero-order, R2)--> [Central (V)] --CL + Vmax/Km-->
;;   (1-DF) x Dose at t=DF x Dose/R2 --> [Depot] --Ka--^
;;
;; Data set / dosing requirements
;;   TWO dose records per administration, both with the full AMT:
;;     CMT=2 (central), RATE=-1  -> DF     x Dose, zero-order at rate R2
;;     CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, released at
;;                                  t = DF x Dose / R2, absorbed with Ka
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM 1CMT seq_fo_zo_rate | Nonlinear and linear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=2
        COMP=(DEPOT, DEFDOSE)
        COMP=(CENTRAL, DEFOBS)

; compartment 1 = depot
; compartment 2 = central (observation)

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V     = THETA(2)*EXP(ETA(2))
  VM    = THETA(3)*EXP(ETA(3))
  KM    = THETA(4)*EXP(ETA(4))
  KA    = THETA(5)*EXP(ETA(5))
  R2    = THETA(6)*EXP(ETA(6))
  LGTDF = LOG(THETA(7)/(1-THETA(7))) + ETA(7)
  DF    = EXP(LGTDF)/(1+EXP(LGTDF))

; ---- structural / input specification ---------------------------------
  F1 = 1 - DF
  F2 = DF
  IF(NEWIND.NE.2) DOSE = 0
  IF(AMT.GT.0)    DOSE = AMT
  ALAG1 = DF*DOSE/R2      ; first-order input starts when the
                          ; zero-order infusion has finished

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V

$DES
  CONC    = A(2)/V
  DADT(1) = -KA*A(1)
  DADT(2) =  KA*A(1) - (CL*CONC + VM*CONC/(KM+CONC))

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- initial estimates -------------------------------------------------
$THETA  (0, 5.0)             ; 1 CL     linear clearance (L/h)
$THETA  (0, 50.0)            ; 2 V      central volume (L)
$THETA  (0, 100.0)           ; 3 VM     Vmax (mg/h)
$THETA  (0, 5.0)             ; 4 KM     Km (mg/L)
$THETA  (0, 1.0)             ; 5 KA     absorption rate constant (1/h)
$THETA  (0, 50.0)            ; 6 R2     zero-order rate into cmt 2 (mg/h)
$THETA  (0.001, 0.5, 0.999)  ; 7 DF     dose fraction (logit-scale IIV)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V
$OMEGA  0.09             ; 3 IIV VM
$OMEGA  0 FIX            ; 4 IIV KM    (free to estimate if supported)
$OMEGA  0.09             ; 5 IIV KA
$OMEGA  0 FIX            ; 6 IIV R2    (free to estimate if supported)
$OMEGA  0 FIX            ; 7 IIV DF    (free to estimate if supported)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=1cmt_seq_fo_zo_rate_mmlin.tab
$TABLE   ID CL V VM KM KA R2 DF ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7
         FIRSTONLY ONEHEADER NOPRINT FILE=1cmt_seq_fo_zo_rate_mmlin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=1cmt_seq_fo_zo_rate_mmlin.tab


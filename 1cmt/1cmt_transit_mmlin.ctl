;; ==========================================================================
;; 1CMT_TRANSIT_MMLIN
;; #21 of the one-compartment set
;;   Transit compartment absorption
;;   Nonlinear and linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --> [T1] --Ktr--> ... --Ktr--> [Tn] --Ktr--> [Depot] --Ka-->
;;     [Central (V)] --CL + Vmax/Km-->
;;           N transit compartments,  Ktr = (N+1)/MTT
;;
;; Data set / dosing requirements
;;   One dose record per administration: CMT=1 (depot), RATE=0.
;;   NOTE: the analytical transit solution tracks only the MOST RECENT
;;   dose.  It is exact for single-dose data; for multiple dosing either
;;   code the transit chain explicitly with a fixed N, or use the
;;   superposition form.
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM 1CMT transit | Nonlinear and linear elimination

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
  MTT   = THETA(6)*EXP(ETA(6))
  NN    = THETA(7)

; ---- structural / input specification ---------------------------------
; --- transit chain (Savic et al. 2007) ------------------------------
;     N transit compartments + depot,  Ktr = (N+1)/MTT
;     the chain is solved analytically, so F1 = 0 and the drug enters
;     the depot through the gamma-density input term coded in $DES.
  KTR   = (NN+1)/MTT
  LNFAC = LOG(2.5066)+(NN+0.5)*LOG(NN)-NN+LOG(1+1/(12*NN))  ; log(N!)
  IF(NEWIND.NE.2) THEN
    TDOS = 0
    DOSE = 0
  ENDIF
  IF(AMT.GT.0) THEN
    TDOS = TIME
    DOSE = AMT
  ENDIF
  F1 = 0

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V

$DES
  TAD = T - TDOS
  INP = 0
  IF(TAD.GT.0.AND.DOSE.GT.0) THEN
    LINP = LOG(DOSE)+LOG(KTR)+NN*LOG(KTR*TAD)-KTR*TAD-LNFAC
    INP  = EXP(LINP)
  ENDIF
  CONC    = A(2)/V
  DADT(1) = INP - KA*A(1)
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
$THETA  (0, 2.0)             ; 6 MTT    mean transit time (h)
$THETA  (0.1, 5.0)           ; 7 NN     number of transit compartments

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V
$OMEGA  0.09             ; 3 IIV VM
$OMEGA  0 FIX            ; 4 IIV KM    (free to estimate if supported)
$OMEGA  0.09             ; 5 IIV KA
$OMEGA  0.09             ; 6 IIV MTT

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=1cmt_transit_mmlin.tab
$TABLE   ID CL V VM KM KA MTT NN ETA1 ETA2 ETA3 ETA4 ETA5 ETA6
         FIRSTONLY ONEHEADER NOPRINT FILE=1cmt_transit_mmlin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=1cmt_transit_mmlin.tab


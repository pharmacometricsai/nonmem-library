;; ==========================================================================
;; 1CMT_TRANSIT_MM
;;   Transit compartment absorption
;;   Nonlinear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --> [T1] --Ktr--> ... --Ktr--> [Tn] --Ktr--> [Depot] --Ka-->
;;     [Central (V)] --Vmax/Km-->
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

$PROBLEM 1CMT transit | Nonlinear elimination

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
  VM    = THETA(1)*EXP(ETA(1))
  KM    = THETA(2)*EXP(ETA(2))
  V     = THETA(3)*EXP(ETA(3))
  KA    = THETA(4)*EXP(ETA(4))
  MTT   = THETA(5)*EXP(ETA(5))
  NN    = THETA(6)

; ---- random effects exported to $TABLE --------------------------------
  ET1   = ETA(1)
  ET2   = ETA(2)
  ET3   = ETA(3)
  ET4   = ETA(4)
  ET5   = ETA(5)

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
  DADT(2) =  KA*A(1) - (VM*CONC/(KM+CONC))

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- residual error variances exported to $TABLE ----------------------
  SG1   = SIGMA(1,1)          ; proportional
  SG2   = SIGMA(2,2)          ; additive

; ---- initial estimates -------------------------------------------------
$THETA  (0, 100.0)           ; 1 VM     Vmax (mg/h)
$THETA  (0, 5.0)             ; 2 KM     Km (mg/L)
$THETA  (0, 50.0)            ; 3 V      central volume (L)
$THETA  (0, 1.0)             ; 4 KA     absorption rate constant (1/h)
$THETA  (0, 2.0)             ; 5 MTT    mean transit time (h)
$THETA  (0.1, 5.0)           ; 6 NN     number of transit compartments

$OMEGA  0.09             ; 1 IIV VM
$OMEGA  0 FIX            ; 2 IIV KM    (free to estimate if supported)
$OMEGA  0.09             ; 3 IIV V
$OMEGA  0.09             ; 4 IIV KA
$OMEGA  0.09             ; 5 IIV MTT

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES PRED
         RES WRES ET1 ET2 ET3 ET4 ET5 SG1 SG2
         ONEHEADER NOPRINT FILE=1cmt_transit_mm.tab
$TABLE   ID VM KM V KA MTT NN ET1 ET2 ET3 ET4 ET5
         FIRSTONLY ONEHEADER NOPRINT FILE=1cmt_transit_mm.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=1cmt_transit_mm_sim.tab


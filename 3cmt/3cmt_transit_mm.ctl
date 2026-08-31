;; ==========================================================================
;; 3CMT_TRANSIT_MM
;; #20 of the three-compartment set
;;   Transit compartment absorption
;;   Nonlinear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --> [T1] --Ktr--> ... --Ktr--> [Tn] --Ktr--> [Depot] --Ka-->
;;        [Central (V2)]     N transit cmts,  Ktr = (N+1)/MTT
;;   [Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
;;   elimination from central: Vmax/Km
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

$PROBLEM 3CMT transit | Nonlinear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=4
        COMP=(DEPOT, DEFDOSE)
        COMP=(CENTRAL, DEFOBS)
        COMP=(PERIPH1)
        COMP=(PERIPH2)

; compartment 1 = depot
; compartment 2 = central (observation)
; compartment 3 = peripheral 1
; compartment 4 = peripheral 2

$PK
; ---- typical values and between-subject variability ------------------
  VM    = THETA(1)*EXP(ETA(1))
  KM    = THETA(2)*EXP(ETA(2))
  V2    = THETA(3)*EXP(ETA(3))
  Q3    = THETA(4)*EXP(ETA(4))
  V3    = THETA(5)*EXP(ETA(5))
  Q4    = THETA(6)*EXP(ETA(6))
  V4    = THETA(7)*EXP(ETA(7))
  KA    = THETA(8)*EXP(ETA(8))
  MTT   = THETA(9)*EXP(ETA(9))
  NN    = THETA(10)

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
  S2 = V2

$DES
  TAD = T - TDOS
  INP = 0
  IF(TAD.GT.0.AND.DOSE.GT.0) THEN
    LINP = LOG(DOSE)+LOG(KTR)+NN*LOG(KTR*TAD)-KTR*TAD-LNFAC
    INP  = EXP(LINP)
  ENDIF
  CONC    = A(2)/V2
  CP1     = A(3)/V3
  CP2     = A(4)/V4
  DIST1   = Q3*(CONC - CP1)
  DIST2   = Q4*(CONC - CP2)
  ELR     = VM*CONC/(KM+CONC)
  DADT(1) = INP - KA*A(1)
  DADT(2) =  KA*A(1) - DIST1 - DIST2 - ELR
  DADT(3) =  DIST1
  DADT(4) =  DIST2

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- initial estimates -------------------------------------------------
$THETA  (0, 100.0)           ; 1 VM     Vmax (mg/h)
$THETA  (0, 5.0)             ; 2 KM     Km (mg/L)
$THETA  (0, 50.0)            ; 3 V2     central volume (L)
$THETA  (0, 10.0)            ; 4 Q3     inter-cmt clearance, periph.1 (L/h)
$THETA  (0, 100.0)           ; 5 V3     peripheral volume 1 (L)
$THETA  (0, 2.0)             ; 6 Q4     inter-cmt clearance, periph.2 (L/h)
$THETA  (0, 200.0)           ; 7 V4     peripheral volume 2 (L)
$THETA  (0, 1.0)             ; 8 KA     absorption rate constant (1/h)
$THETA  (0, 2.0)             ; 9 MTT    mean transit time (h)
$THETA  (0.1, 5.0)           ; 10 NN     number of transit compartments

$OMEGA  0.09             ; 1 IIV VM
$OMEGA  0 FIX            ; 2 IIV KM    (free it if the data support it)
$OMEGA  0.09             ; 3 IIV V2
$OMEGA  0.09             ; 4 IIV Q3
$OMEGA  0.09             ; 5 IIV V3
$OMEGA  0.09             ; 6 IIV Q4
$OMEGA  0.09             ; 7 IIV V4
$OMEGA  0.09             ; 8 IIV KA
$OMEGA  0.09             ; 9 IIV MTT

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=3cmt_transit_mm.tab
$TABLE   ID VM KM V2 Q3 V3 Q4 V4 KA MTT NN ETA1 ETA2 ETA3 ETA4 ETA5 ETA6
         ETA7 ETA8 ETA9
         FIRSTONLY ONEHEADER NOPRINT FILE=3cmt_transit_mm.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=3cmt_transit_mm.tab


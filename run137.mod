;; 1. Based on: run135
;; 2. Description: IV+IG, Three dose, CL, BIO No nonlinear PK
;; x1. Author: RZ

$PROBLEM  METABOLISM OF ARCTIGENIN IN RATS AFTER IV AND IG - ADVAN6
$DATA   Arctigenin_IG+IV_for_nonmem.csv IGNORE=@
$INPUT  ID TIME DV DOSE=AMT MDV CMT EVID DL RTE
$SUBROUTINE ADVAN6 TRANS1 TOL=6
$MODEL
COMP = (LUMEN DEFDOSE)			; 1. intestine lumen
COMP = (AR_CENTR)			; 2. AR in central
COMP = (AG_CENTR)			; 3. AG in central
COMP = (AA_CENTR)			; 4. AA in Central
;COMP = (AG_INTES)			; 5. AG from Intestine

$PK
KA = THETA(1)*EXP(ETA(1))		; absorption rate of AR
K20= THETA(2)*EXP(ETA(2))		; elimination rate of AG from Central
K30= THETA(3)*EXP(ETA(3))		; elimination rate of AA from central comp
V2 = THETA(4)*EXP(ETA(4))		; volume of central
S2 = V2
S3 = V2
S4 = V2
S5 = V2
CLI= THETA(5)*EXP(ETA(5))		; Clearance of intestine
CLB= THETA(6)*EXP(ETA(6))		; Clearance of Blood
CLL= THETA(7)*EXP(ETA(7))		; Clearance of Liver
F1 = THETA(8)*EXP(ETA(8))		; Bioavailability
tKi = CLI/V2
Kb = CLB/V2
Kl = CLL/V2

$DES
;tKi=Ki
DADT(1)=-KA*A(1)
IF(RTE==1)THEN	; IV
DADT(2) = KA*A(1)-Kl*A(2)-Kb*A(2)
DADT(3) = Kl*A(2)-K20*A(3)
ELSE			; IG
DADT(2) = KA*A(1)-Kb*A(2)-tKi*A(2)
DADT(3) = tKi*A(2)-K20*A(3)	
ENDIF
DADT(4) = Kb*A(2)-K30*A(4)

$ERROR
Y = F+ERR(1)

$THETA
(0, 9) ; KA(1/h)
(0, 0.45, ) ; K20(1/h)
(0, 0.8, 3) ; K30(1/h)
(0, 1) ; V2(L)
(0, 15) ; CLI(L/h)
(0, 10) ; CLB(L/h)
(0, 10) ; CLL(L/h)
(0,0.01,1) ; F1

$OMEGA
0 FIX ; IIV KA
0 FIX ; IIV K20
0 FIX ; IIV K30
0 FIX ; IIV V2
0 FIX ; IIV CLI
0 FIX ; IIV CLB
0 FIX ; IIV CLL
0 FIX ; IIV BIO

$SIGMA
0.000432 ; Residual error

$EST METHOD=0 MAXEVAL=999 NOABORT SIG=5 PRINT=1
$COV

; Xpose
$TABLE ID TIME DV MDV CMT DL ONEHEADER NOPRINT FILE=sdtab137.tab
$TABLE ID Ka K20 K30 V2 ONEHEADER NOPRINT FILE=patab137.tab


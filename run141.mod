;; 1. Based on: run140
;; 2. Description: IG, Three dose, CL, BIO, CL adjust by dose
;; x1. Author: RZ

$PROBLEM  METABOLISM OF ARCTIGENIN IN RATS AFTER IV - ADVAN6
$DATA   Arctigenin_IG_for_nonmem.csv IGNORE=@
$INPUT  ID TIME DV DOSE=AMT MDV CMT DL
$SUBROUTINE ADVAN6 TOL=6
$MODEL
COMP=(LUMEN DEFDOSE)			; 1. intestine lumen
COMP=(AG_CENTR)				; 2. AG in central
COMP=(AA_CENTR)				; 3. AA in Central
COMP=(AR_INTES)				; 4. AR in Intestine compartment
COMP=(AG_INTES)				; 5. AG in intestine
COMP=(AA_INTES)				; 6. AA in intestine


$PK
Ka = THETA(1)*EXP(ETA(1))		; absorption rate of AR
K20= THETA(2)*EXP(ETA(2))		; elimination rate of AG from Central
K30= THETA(3)*EXP(ETA(3))		; elimination rate of AA from central comp
V2=  THETA(4)				; volume of central
S2=  V2
S3=  V2
CLI= THETA(5)*EXP(ETA(4))		; Clearance of intestine
CLB= THETA(6)*EXP(ETA(5))		; Clearance of Blood
F1= THETA(7)*EXP(ETA(6))		; Bioavailability
Ki=CLI/V2
Kb=CLB/V2
;K52= THETA(8)
;K63= THETA(9)
;K63= THETA(10)
;K36= THETA(11)


$DES
IF(DL==1)THEN
tKi=Ki*1.314277
ENDIF
IF(DL==3) THEN
tKi=Ki*0.610052
ENDIF
IF(DL==2)THEN
tKi=Ki
ENDIF

DADT(1)=-Ka*A(1)
DADT(2)= tKi*A(5)-K20*A(2)
DADT(3)= Kb*A(6)-K30*A(3)
DADT(4)= Ka*A(1)-tKi*A(4)-Kb*A(4)
DADT(5)= tKi*A(4)-tKi*A(5)
DADT(6)= Kb*A(4)-Kb*A(6)

$ERROR
IPRED = F
;W =SQRT(THETA(8)**2*IPRED**2 + THETA(9)**2)
W =1
Y = IPRED + W*EPS(1)
IRES = DV-IPRED
IWRES = IRES/W

$THETA
(0, 40.6) ; Ka
(0, 0.45) ; K20
(0, 0.849) ; K30
(0, 0.67) ; V2(L)
(0, 9.87) ; CLI
(0, 11.6) ; CLB
(0, 0.01,1) ; BIO

$OMEGA
 0 FIX  ; IIV Ka
 0 FIX  ; IIV K20
 0 FIX  ; IIV K30
 0 FIX  ; IIV CLI
 0 FIX  ; IIV CLB
 0 FIX  ; IIV BIO

$SIGMA
 0.000172 ; Residual error

$EST METHOD=0 MAXEVAL=999 NOABORT SIG=5 PRINT=1
$COV

; Xpose
$TABLE ID TIME DV MDV CMT DL ONEHEADER NOPRINT FILE=sdtab141.tab
$TABLE ID Ka K20 K30 V2 ONEHEADER NOPRINT FILE=patab141.tab


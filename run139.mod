;; 1. Based on: run134
;; 2. Description: IG, Three dose, CL, BIO
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
;IF(DL==1)THEN
;tKi=Ki*1.314277
;tKb=Kb*1.209045
;ENDIF
;IF(DL==2)THEN
tKi=Ki
;ENDIF
;IF(DL==3)THEN
;tKi=Ki*0.610052
;tKb=Kb*0.902952
;ENDIF
tKb=Kb

DADT(1)=-Ka*A(1)
DADT(2)= tKi*A(5)-K20*A(2)
DADT(3)= tKb*A(6)-K30*A(3)
DADT(4)= Ka*A(1)-tKi*A(4)-tKb*A(4)
DADT(5)= tKi*A(4)-tKi*A(5)
DADT(6)= tKb*A(4)-tKb*A(6)

$ERROR
Y = F+ERR(1)

$THETA
(0, 35) ; Ka
(0, 0.45, ) ; K20
(0, 0.8, 3) ; K30
(0, 1) ; V2
(0, 15) ; CLI
(0, 10) ; CLB
(0,0.01,1) ; BIO
;(0, 9) ; K52
;(0, 9) ; K25
;(0, 9) ; K63
;(0, 9) ; K36

$OMEGA
0 FIX ; IIV Ka
0 FIX ; IIV K20
0 FIX ; IIV K30
0 FIX ; IIV CLI
0 FIX ; IIV CLB
0 FIX ; IIV BIO

$SIGMA
0.000432 ; Residual error

$EST METHOD=0 MAXEVAL=999 NOABORT SIG=5 PRINT=1
$COV

; Xpose
$TABLE ID TIME DV MDV CMT DL ONEHEADER NOPRINT FILE=sdtab139.tab
$TABLE ID Ka K20 K30 V2 ONEHEADER NOPRINT FILE=patab139.tab


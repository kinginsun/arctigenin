;; 1. Based on: run144
;; 2. Description: FINAL IG, Three dose, CL, BIO, Constant Km
;; x1. Author: RZ

$PROBLEM  METABOLISM OF ARCTIGENIN IN RATS AFTER IV - ADVAN6
$DATA   Arctigenin_IG_for_nonmem.csv IGNORE=@
$INPUT  ID TIME DV DOSE=AMT MDV CMT DL
$SUBROUTINE ADVAN6 TOL=6
$MODEL
COMP=(LUMEN DEFDOSE)			; 1. intestine lumen
COMP=(AG_CENTR)				; 2. AG in central
COMP=(AA_CENTR)				; 3. AA in Central
COMP=(INTESTIN)				; 4. Intestine compartment
;COMP=(AG_INTES)			; 5. AG in intestine
;COMP=(AG_PER)				; 5. AG in perp
;COMP=(AA_PER)				; 6. AA in perp

$PK
Ka = THETA(1)*EXP(ETA(1))		; absorption rate of AR
K2e= THETA(2)*EXP(ETA(2))		; elimination rate of AG from Central
K3e= THETA(3)*EXP(ETA(3))		; elimination rate of AA from central comp
V2=  THETA(4)				; volume of central
S2=  V2
S3=  V2
CLI= THETA(5)*EXP(ETA(4))		; Clearance of intestine
CLB= THETA(6)*EXP(ETA(5))		; Clearance of Blood
F1= THETA(7)*EXP(ETA(6))		; Bioavailability
Ki=CLI
Kb=CLB
K20=K2e/1000
K30=K3e/1000
tKi=Ki

$DES
;IF(DL==1) tKi=Ki*1.314277
;IF(DL==2) tKi=Ki
;IF(DL==3) tKi=Ki*0.610052

DADT(1)=-Ka*A(1)
DADT(2)= tKi*A(4)-K20*A(2)
DADT(3)= Kb*A(4)-K30*A(3)
DADT(4)= Ka*A(1)-tKi*A(4)-Kb*A(4)

$ERROR
Y = F+ERR(1)

$THETA
(0, 0.3) ; Ka(1/min)
(0, 1) ; K2e(1/min)
(0, 10) ; K3e(1/min)
(0, 0.289 FIX) ; V2(L)
(0, 0.7) ; Ki (1/min)
(0, 0.8) ; Kb (1/min)
(0,0.01,1) ; BIO

$OMEGA
0 FIX ; IIV Ka
0 FIX ; IIV K20
0 FIX ; IIV K30
0 FIX ; IIV Ki
0 FIX ; IIV Kb
0 FIX ; IIV BIO

$SIGMA
0.000432 ; Residual error

$EST METHOD=0 MAXEVAL=999 NOABORT SIG=5 PRINT=1
$COV

; Xpose
$TABLE ID TIME DV MDV CMT DL ONEHEADER NOPRINT FILE=sdtab145.tab
$TABLE ID Ka K20 K30 V2 ONEHEADER NOPRINT FILE=patab145.tab

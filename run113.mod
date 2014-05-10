;; 1. Based on: run101
;; 2. Description: FINAL IV Arctigenin to AA + AG, IV, Adjust Vm and Km
;; x1. Author: RZ

$PROBLEM  METABOLISM OF ARCTIGENIN IN RATS AFTER IV - ADVAN6
$DATA   Arctigenin_IV_for_nonmem.csv IGNORE=@
$INPUT  ID TIME DV DOSE=AMT MDV CMT DL
$SUBROUTINE ADVAN6 TOL=6
$MODEL
COMP=(CENTRAL,DEFDOSE,DEFOBS)		; 1. Arctigenin (AR) central
COMP=(AG_LIVER)				; 2. AG Liver, equal to central
COMP=(AA_CENTR)				; 3. AA Central

$PK
K10= THETA(1)*EXP(ETA(1))		; elimination rate of AR from central comp
K20= THETA(2)*EXP(ETA(2))		; elimination rate of AG from Liver, not considering the transfer from Liver to Central
K30= THETA(3)*EXP(ETA(3))		; elimination rate of AA from central comp
V1=  THETA(4)*EXP(ETA(4))
S1=  V1
S2=  V1
S3=  V1
Fi1= THETA(5)*EXP(ETA(5))		; scale factor for blood hydrolysis
Fi2= THETA(6)*EXP(ETA(6))		; scale factor for liver glucuronidation
Km1= THETA(7)				; AA blood hydrolysis
Vm1= THETA(8)
Km2= THETA(9)				; AG liver glucuronidation
Vm2= THETA(10)
$DES
DADT(1)=-K10*A(1)-Fi1*Vm1*A(1)/(Km1*V1+A(1))-Fi2*Vm2*A(1)/(Km2*V1+A(1))
DADT(2)= Fi2*Vm2*A(1)/(Km2*V1+A(1))-K20*A(2)
DADT(3)= Fi1*Vm1*A(1)/(Km1*V1+A(1))-K30*A(3)

$ERROR
IPRED = F
W =1; SQRT(THETA(11)**2*IPRED**2 + THETA(12)**2)
Y = IPRED + W*EPS(1)
IRES = DV-IPRED
IWRES = IRES/W

$THETA
(0 FIX) ; K10
(0, 0.206) ; K20
(0, 0.132) ; K30
(0, 0.311) ; V1
(0, 1.6) ; Fi1
(0, 0.0638) ; Fi2
(13.53) FIX ; Km1
(0.292) FIX ; Vm1
(230) FIX ; Km2
(1.23) FIX ; Vm2
;(0, 0.1)  ; Prop.RE (sd)
;(0, 0.102) ; Add.RE (sd)

$OMEGA
0 FIX ; IIV K10
0 FIX ; IIV K20
0 FIX ; IIV K30
0 FIX ; IIV V1
;(0.1) ; IIV V2
;(0.1) ; IIV V3
0 FIX ; IIV Fi1
0 FIX ; IIV Fi2

$SIGMA
0.326 ; Residual error

$EST METHOD=0 MAXEVAL=999 NOABORT SIG=5 PRINT=1 POSTHOC
$COV

; Xpose
$TABLE ID TIME DV MDV CMT DL IPRED IWRES ONEHEADER NOPRINT FILE=sdtab113.tab
$TABLE ID K10 K20 K30 V1 Fi1 Fi2 ONEHEADER NOPRINT FILE=patab113.tab


;; 1. Based on: run129
;; 2. Description: FINAL IV Arctigenin to AA + AG, IV, CLL CLB
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
K10= THETA(1)*EXP(ETA(1))		; elimination rate of AR from central comp (1/min)
K20= THETA(2)*EXP(ETA(2))		; elimination rate of AG from Liver, not considering the transfer from Liver to Central
K30= THETA(3)*EXP(ETA(3))		; elimination rate of AA from central comp
V1=  THETA(4)*EXP(ETA(4))		; unit: mL
S1=  V1/1000
S2=  V1/1000
S3=  V1/1000
CLB=  THETA(5)*EXP(ETA(5))		; clearance of Blood (mL/min)
CLL=  THETA(6)*EXP(ETA(6))		; clearance of Liver (mL/min)
KB=CLB/V1
KL=CLL/V1

$DES

DADT(1)=-K10*A(1)-KB*A(1)-KL*A(1)
DADT(2)= KL*A(1)-K20*A(2)
DADT(3)= KB*A(1)-K30*A(3)

$ERROR
Y = F + EPS(1)

$THETA
(0 FIX)    ; K10
(0, 0.206) ; K20
(0, 0.132) ; K30
(0, 300)   ; V1(mL)
(0, 25)   ; CLB(mL/min)
(0, 10)   ; CLL(mL/min)

$OMEGA
0 FIX ; IIV K10
0 FIX ; IIV K20
0 FIX ; IIV K30
0 FIX ; IIV V1
0 FIX ; IIV CLB
0 FIX ; IIV CLL

$SIGMA
0.01 ; Residual error

$EST METHOD=0 MAXEVAL=999 NOABORT SIG=5 PRINT=1 POSTHOC
$COV

; Xpose
$TABLE ID TIME DV MDV CMT DL ONEHEADER NOPRINT FILE=sdtab130.tab


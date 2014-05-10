;; 1. Based on: run118
;; 2. Description: FINAL IG Arctigenin to AA + AG, IG, mean Cl 2, All, Adding BIO
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
COMP=(AG_PER)				; 5. AG in perp
COMP=(AA_PER)				; 6. AA in perp

$PK
Ka = THETA(1)*EXP(ETA(1))		; absorption rate of AR
K20= THETA(2)*EXP(ETA(2))		; elimination rate of AG from Central
K30= THETA(3)*EXP(ETA(3))		; elimination rate of AA from central comp
;K52= THETA(4)				; AG from intestine to central
;K63= THETA(5)				; AR from intestine to central
V2=  THETA(4)				; volume of central
;V4=  THETA(7)				; volume of intestine
S2=  V2
S3=  V2
Ki= THETA(5)				; scale factor for blood hydrolysis
Kb= THETA(6)				; scale factor for liver glucuronidation
;Fi3= THETA(10)				; scale factor for intestine glucuronidation
;Km1= THETA(11)				; AA blood hydrolysis
;Vm1= THETA(12)
;Km2= THETA(13)				; AG liver glucuronidation
;Vm2= THETA(14)
;Km3= THETA(15)				; AG intestine glucuronidation
;Vm3= THETA(16)
K52= THETA(7)
K25= THETA(8)
K63= THETA(9)
K36= THETA(10)
BIO= THETA(11)
$DES
DADT(1)=-Ka*A(1)
DADT(2)= Ki*A(4)-K20*A(2)-K25*A(2)+K52*A(5)
DADT(3)= Kb*A(4)-K30*A(3)-K36*A(3)+K63*A(6)
DADT(4)= BIO*Ka*A(1)-Ki*A(4)-Kb*A(4)
DADT(5)= K25*A(2)-K52*A(5)
DADT(6)= K36*A(3)-K63*A(6)

$ERROR
Y = F+ERR(1)

$THETA
(0, 9) ; Ka
(0, 0.45) ; K20
(0, 2) ; K30
(0, 16) ; V2
(0, 1.5) ; Ki
(0, 10) ; Kb
(0, 1) ; K52
(0, 0.225) ; K25
(0, 1.36) ; K63
(0, 7.45) ; K36
(0,0.01,0.1) ; BIO
$OMEGA
0 FIX 
0 FIX
0 FIX

$SIGMA
0.000432 ; Residual error

$EST METHOD=0 MAXEVAL=999 NOABORT SIG=5 PRINT=1
$COV

; Xpose
$TABLE ID TIME DV MDV CMT DL ONEHEADER NOPRINT FILE=sdtab128.tab
$TABLE ID Ka K20 K30 V2 ONEHEADER NOPRINT FILE=patab128.tab


;; 1. Based on: run138
;; 2. Description: IG, Three dose, CL, BIO, 7 CMTs
;; x1. Author: RZ

$PROBLEM  METABOLISM OF ARCTIGENIN IN RATS AFTER IV AND IG - ADVAN6
$DATA   Arctigenin_IG_for_nonmem.csv IGNORE=@
$INPUT  ID TIME DV DOSE=AMT MDV CMT DL
$SUBROUTINE ADVAN6 TRANS1 TOL=6
$MODEL
COMP = (LUMEN DEFDOSE)			; 1. intestine lumen
COMP = (AR_CENTR)			; 2. AR in central
COMP = (AG_CENTR)			; 3. AG in central
COMP = (AA_CENTR)			; 4. AA in Central
COMP = (AR_INTES)			; 5. AR in Intestine
COMP = (AG_INTES)			; 5. AG in Intestine
COMP = (AA_INTES)			; 5. AA in Intestine

$PK
KA = THETA(1)*EXP(ETA(1))		; absorption rate of AR

K20= THETA(2)*EXP(ETA(2))		; elimination rate of AG from Central
K30= THETA(3)*EXP(ETA(3))		; elimination rate of AA from central comp
V2 = THETA(4)*EXP(ETA(4))		; volume of central
V5 = THETA(11)*EXP(ETA(11))		; volume of intestine
S2 = V2
S3 = V2
S4 = V2
;S5 = V2
Ki= THETA(5)*EXP(ETA(5))		; Clearance of intestine
Kb= THETA(6)*EXP(ETA(6))		; Clearance of Blood
Kl= THETA(7)*EXP(ETA(7))		; Clearance of Liver
F1 = THETA(8)*EXP(ETA(8))		; Bioavailability
K10= THETA(9)*EXP(ETA(9))		; elimination rate of AR from Central
K52= THETA(10)*EXP(ETA(10))		; transfter rate from intestine to Central
K63= THETA(12)*EXP(ETA(12))	
K36= THETA(13)*EXP(ETA(13))	
K74= THETA(14)*EXP(ETA(14))	
K47= THETA(15)*EXP(ETA(15))	

$DES

DADT(1)=-KA*A(1)
DADT(2) = K52*A(5)-Kb*A(2)-Kl*A(2)-K10*A(2)
DADT(3) = K63*A(6)-K36*A(3)-K20*A(3)
DADT(4) = K74*A(7)-K47*A(4)-K30*A(4)
DADT(5) = KA*A(1)-Ki*A(5)-Kb*A(5)-K52*A(5)
DADT(6) = Ki*A(5)-K63*A(6)+K36*A(3)
DADT(7) = Kb*A(5)-K74*A(7)+K47*A(4)

$ERROR
Y = F+ERR(1)

$THETA
(0, 0.3) ; KA(1/min)
(0.18 FIX) ; K20(1/min)
(0.14 FIX) ; K30(1/min)
(0.289 FIX) ; V2(L)
(0, 1) ; Ki(1/min)
(0, 0.101 FIX) ; Kb(1/min)
(0, 0.0244 FIX) ; Kl(1/min)
(0,0.01,1) FIX ; F1
0 FIX  ; K10
0 FIX  ; K52
(0,1) ; V5(L)
(0,0.1) ; K63
(0,0.1) ; K36
(0,0.1) ; K74
(0,0.1) ; K47

$OMEGA
0 FIX ; IIV KA
0 FIX ; IIV K20
0 FIX ; IIV K30
0 FIX ; IIV V2
0 FIX ; IIV CLI
0 FIX ; IIV CLB
0 FIX ; IIV CLL
0 FIX ; IIV BIO
0 FIX ; IIV K10
0 FIX ; IIV K52
0 FIX ; IIV V5
0 FIX ; IIV K63
0 FIX ; IIV K36
0 FIX ; IIV K74
0 FIX ; IIV K47

$SIGMA
0.000432 ; Residual error

$EST METHOD=0 MAXEVAL=999 NOABORT SIG=5 PRINT=1
$COV

; Xpose
$TABLE ID TIME DV MDV CMT DL ONEHEADER NOPRINT FILE=sdtab142.tab
$TABLE ID Ka K20 K30 V2 ONEHEADER NOPRINT FILE=patab142.tab


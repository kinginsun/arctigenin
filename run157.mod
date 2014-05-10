;; 1. Based on: run156
;; 2. Description: Good IV+IG, Three dose, integrating MM functions
;; x1. Author: user

$PROBLEM  METABOLISM OF ARCTIGENIN IN RATS AFTER IV AND IG - ADVAN6
$DATA   Arctigenin_IG+IV_for_nonmem.csv IGNORE=@
$INPUT  ID TIME DV DOSE=AMT MDV CMT EVID DL RTE
$SUBROUTINE ADVAN6 TRANS1 TOL=9
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
S2 = V2
S3 = V2
S4 = V2
;S5 = V2
SF1= THETA(5)*EXP(ETA(5))		; Scale factor for intestine metabolism
SF2= THETA(6)*EXP(ETA(6))		; Scale factor for blood hydrolysis
SF3= THETA(7)*EXP(ETA(7))		; Scale factor for liver metabolism
F1 = THETA(8)*EXP(ETA(8))		; Bioavailability
K10= THETA(9)*EXP(ETA(9))		; elimination rate of AR from Central
K52= THETA(10)*EXP(ETA(10))	

; Up-Scaled Km and Vm from in -vitro experiments
Km1= THETA(11)				; AG intestine glucuronidation 
Vm1= THETA(12)
Km2= THETA(13)				; AA blood hydrolysis 
Vm2= THETA(14)
Km3= THETA(15)				; AG liver glucuronidation
Vm3= THETA(16)
;V5 = THETA(17)
$DES

Ki=SF1*Vm1/(Km1*V2+A(5))
Kb=SF2*Vm2/(Km2*V2+A(2))
Kl=SF3*Vm3/(Km3*V2+A(2))
DADT(1)=-KA*A(1)
DADT(2) = K52*A(5)-Kb*A(2)-Kl*A(2)-K10*A(2)
DADT(3) = Ki*A(6)+Kl*A(2)-K20*A(3)
DADT(4) = Kb*A(7)+Kb*A(2)-K30*A(4)
DADT(5) = KA*A(1)-Ki*A(5)-Kb*A(5)-K52*A(5)
DADT(6) = Ki*A(5)-Ki*A(6)
DADT(7) = Kb*A(5)-Kb*A(7)

$ERROR
Y = F+ERR(1)

$THETA
(0, 0.03) ; KA(1/min)
(0, 0.1) ; K20(1/min)
(0, 0.1) ; K30(1/min)
(0, 0.3) ; V2(L)
(0, 1) ; SF1
(0, 1) ; SF2
(0, 1) ; SF3
(0, 0.01,1) ; F1
(0 FIX)  ; K10
(0 FIX) ; K52
(203.8) FIX ; Km1(umol/L)
(12.449) FIX ; Vm1(umol/min)
(13.53) FIX ; Km2
(0.292) FIX ; Vm2
(230.2) FIX ; Km3
(1.231) FIX ; Vm3
;(0, 0.1)    ; V5(L)

$OMEGA
0 FIX  ; IIV KA
0 FIX  ; IIV K20
0 FIX  ; IIV K30
0 FIX  ; IIV V2
0 FIX  ; IIV Ki
0 FIX  ; IIV Kb
0 FIX  ; IIV Kl
0 FIX  ; IIV BIO
0 FIX  ; IIV K10
0 FIX  ; IIV K52

$SIGMA
0.0186 ; Residual error

$EST METHOD=0 MAXEVAL=999 NOABORT SIG=5 PRINT=1
$COV

; Xpose
$TABLE ID TIME DV MDV CMT DL RTE ONEHEADER NOPRINT FILE=sdtab157.tab
$TABLE ID Ka K20 K30 V2 ONEHEADER NOPRINT FILE=patab157.tab


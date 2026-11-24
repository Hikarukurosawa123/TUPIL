#!/usr/bin/env python3
"""Run kidney-ROI H-scan extraction as one self-contained Python script.

Defaults reproduce the existing sigma_auto1e6 dataset: TGC removal, kidney
attenuation 0.5, non-kidney attenuation 0.8, attenuation taper enabled, and
probe-derived sigma (1.0e6 Hz for C3).

This file embeds config.py, core.py, hscan.py, kidney_aec_202603.py, and
extract_all_hscan_regional_attenuation.py. External requirements are numpy,
scipy, pandas, PyYAML, and matplotlib.
"""

from __future__ import annotations

import argparse
import base64
import os
from pathlib import Path
import sys
import types
import zlib


_BUNDLED_MODULES = {
    "config": (
        "c-pO2ZExE)5dQ98aZoV8PGcu_>U1>_pj&_>1DdQ@g7w=BS~^-"
        "pC{iV<1oE%nk<^RjrbD|;A8dlWd+xbcGD(svt+?1$hB=`XSaeRPtp%<Pe5q{7pk(P2+u|%qlF>*Mr8W*!Rg^U_3#uG(qYKEGW2t1"
        ";qJ>7L+Zv5N^y}~o(dcRf^d*NJqKGGOJAFu*g3{;|ttDQ8gjYI(8&qhR>==h_ZdidU5X!+NOy^#Xz+Bw#CD!bjbho6of{oUags+*"
        "jQ9=s!O@RPN=g^%8p5w*w{l^KUsGO%{cmPteWfF)`%r<C1^MGK%=>`Q!!UJ^HTebJg&Yl-mOx8OD!9nSsz$nJ!b1u>xJ{i4<UWC;"
        "OtTP0ZT%9Lylzf~olSUjq<+T+lku8~zh(i-2XtiAdmzo9M%eQA<*p`TyXKr<+GGYhCscn&Bk`7?*TOSt7i=p$m^(o~-"
        "cCzG9U<;R#X)Fu>ZC=(q?dhD&0+xM2EMbF={uNAv&`&->&5DfpOprNl@Jod%tvf_3wguZqJjaOXn97*RAj4-7%g`tGhE-"
        "H|2Yf>r?a^AymJ#Ht$$!CInupp5T0IRPx!YX;K{;^Sr{*90e)qroeevJ@zW7nUW5H1Eb+r<V{o)*|DjL7i#$f7~L72%hAJ7gR!S)"
        "nVjhu^A5Cv}~`}K5ucf0#*0ckOTBM<*1GKDm(9F)H1wba{C_zsJbddO7<PoSdSp3lx&{ZrEJ39KW;S3DZzKYsu1yD9B~Dy?1t{yB"
        "7uR?~kAuE;5i+k^g-"
        "sJ1TmC4!9R1PK%#!N=#}@%k}n=DKXP=y%)z!7;H~EGHg?B&d>T1`YHI&RsX@&K7Tn{eFG@2&JY|+iSI><Br*<;Q$wQ?Q%3UfY${b"
        "cGEXZ&_`2+@fRdW6c~7+FdWNL*4=FD?u*0Qp{d8UEwpy|>+}Cbs5OPF^=KJVO`Dm1X~XCp=^%m<Szl|Dby4>^rIxkQg`iraBi|2S"
        "sd4EBcSO@Tt}YSA^A$X>XW@2F_4GzdR-eFP1$UTHW1he{E%meavAIPit@}RW1HTm|XkApjt;H}!d+G7Im%Rcd|DM(@@}I(1Op4v`"
        "t*}?qn?h%}U(WCM4#VUhG*OF7ZsHKm(oUC~kas_<yqG~t#G2dZ13e9>OkEsw8FJvo*rMc<ewh0A1(RtdnZv^YDKNU-"
        "jJfW4;Xi{k8Wi=@jBa3h2&Fh4&!--I76$I*IIMYkLaBQ5>V8<<Pey+O@#*Tm"
    ),
    "core": (
        "c-rlKU323|a^U+!%zr4e?I56rAd<a39@9&&`%67KPk%(KXE$uMi~xy(gl!Pu0H7qctcZKLd%cJYKRNE<{=&W8vwzC|f%^-"
        "Wne`2Tlr-HtyLX<BnQjuO%F4>j%CE{Q4hDl~WyM*#E*4Ez7Bx%DimmgeitBP+B&^QPia2Mhs$B58&Wf{eFc^IO)z@ET%T-"
        "x5tgQ9Nihuo8T9r$-iktI1n=|?REqv9DiuH1}VR6lhm26!vvehP(<H%O$Sw81g6R_iR`1(Fu^25XT0XySObeR=Nc^R-vo-"
        "fV=R>v1-Wgewj-ta15OP=6I2XvY?3NCyM@UOJuKhlT#Ji|v>s{W=2PWP--z4J0omT~KC0pqU9JZ`vZCLWCwJ?Bll%#mbYf0ghQ$P"
        "_10m3n|DFXAQlr(b`?;1AIC9lT=i4p|a6F{44|SrJzohHXR9O#J<%I;&+r{G(touWl-TN+f4ZX<(WD2F5$ytgFHX)2vrHpEB<p@Z"
        "c3ovIS730qWOHyj;~RuBvzw2$1yI*QkpR=SdZ1MOv~M^M+Y64A?NO0IwQ9<yj#<)Br4ZK4zcz(AOaBG3978c6}4D*JZ)Q-"
        "!?5SGnj0-"
        ";)NzD0IudkAH}SSFQfUU;r28iFe3D?7`<R}wiY#E3D7LELeI$RF0nkp$ESc9@H<3iiE%tVVu+mfe&q)=PFcYqGay3K<a)AZ5*A)I"
        "!C?jB{4Exlmor#H;U`{|wTI{^N}mGBn$3#OhJe^Nzp{s_KS3vID|m6r_QY4kk@f@R%c<!vS)qsZY*;}Hl>Dx+@@=%n-}F-"
        "eDNh}c^i^D(ad&Ediw%-u-%3uFsiBYC(2(tfVBOBQZ*Sw9H9n6Q=K^o{HPY*38Zh`zj|#-"
        "YLFel<<(22Q$3HDJ>Q8`2x0ORm_q7il>+^WUy^~$+IrZ%|WH?)2#;f=W=z)&{HVW9-"
        "cNw*+U|265eh=!9C+r#<hBV^I<n#t7AZAPaKoJU@HX~?r<vba~6Hg)ruLxmZk-"
        "8!+)vp#>N|hHETeT{@<2fjaFQ;2u5P#9Izi8M$7Y(}xyuSTK!+r@3D@gfS8g{pCwb$97mu5x1nu50-fj^(|qKVQn2TpnD-"
        "Um!Q%#wgzWOdd=NHFy|ODqK}J{<G7TAXtcVZ^FG6F^DDfJ8RuoYgD7$kGf5kYp8KG$knM6)D~GxMARc6BKL-ng}8iwkSar7Z|?CF"
        "^6xj4SF}xwZZQu8eJ<4Z=$_n*sJnsG^y_2!#Hg`N><ZkXex+^JDJddB+(#vjPzb3!LEz!$2Ba(1d)dcS8aF}ATtFUG5QX&KyR7H!"
        "?O@H^P9KNQS%?cKYj<qJu#t9;Wl$70inGRfYNoIYa(I)G^=XKG-Oqa>lsKGhM55?p>cJY)trE-"
        "0b_Vzk0%ED4wN<u^dT6d*JX2v5MET3m=Uiaj?6ct*l9hXj%S=+gN)bl8OLE8hpkBo;EJ+g{3@%}bOE$8_`aB-"
        "z??jS<r?M&%LZaRGLgf#KzC@6mzN;D^Nk}LN6+QrU}B*lg5vwpbLOGC-(STEtoxQas`)CeAiN}iFjHRD3kcb30bv<$*m-"
        "=x8FbBX2|y6NIEM}kpae`kg?|W*JO*4%_t2i#n*s=?sF87WLpZLg8U&)^$={N9Bc4K-"
        "4!skn`%wzCGdxzjtrc&DqN#1pYOU~qJp(Na!E?pW_!SZhIVj!%n}9wA5V%-"
        "yPl<29NdBEQhJ8TA!n79Uas|kGRXX_R9a_Lr0^!4fRERvU`Del%^#V2;zNt3TmNAypDp*;EIY=MyfmCeU5(71$1-"
        "=Qh8XKu^zb~^wDrEu8_enGWYoVq*ITQVU+G7-jG1v}Mm||+wv#zaGVH-"
        "4lwcx9!yHIT^3!66(&s*wth?G9R8MYUwO<0L}+h?^NTd335Y^S_BrY~E9zPO*^3z+ECV*(EWq_AAnM2##A;#mepOlu3R9m1CNv<#"
        "(cYv?5))$6&0tmkOPJbEr6wv}0ZUas>5mwgk%0K+4i#kHCYFu_~IHQ$E?$w4zTfUbYL`FvNFHi8m**A&!uF`92>7h|V(TFBpOZTk"
        "~VqX5u5nj%Z92NA8?;}*+#*tVLQZ?>ESs)G>VMu_WHWAy`=PN`KJDeX4DYnvB6;&CU-"
        ";Ix&K+@T(9g!Q6O+=fv=)qaM<i3``(C~w=#W6G|r!EOeQ@BP41_Z=PIUJC=EQY?#QXgu8Ab;!VpB{vuS22BDpv^|Pz0o4m(0m)op"
        "hj6bw#$WzRxR<q`?ZP=OHbM#TB0+u$Z>Ebn2;i^{c&O-"
        "pciW~|aBmZ*I(97`&Eqo+T?g`?&faMOy}Q}_W0|~&U#;?Nku|L~NQAntnNzSpqX=8qZ^MwGEnW0SvM;D{Oj56f2D0|gSl?@6W5gA"
        "?I<&`$vlrceK%LVj()|=qyZ0B?K3l>6iBwP<pPyF`bu|%&UE-"
        "r44txi_)Sweu7PDcMhE;qyw6qZFA8O`^mPsq=&L|~f#p(I)+6mA|D6Mhi_N685ZX6gyRi8t`OFwmTNw@DG01UGOVaSLu`;<jXSM;"
        "nF?CGu*J<Fa&Lc1J5)EX)3xf8srxVKp8e_7$x&_6e5G_qK%mvb<Ml*k791Hrx8lC~FD4Q5;9@)V8529s#KOQLEmNyTd2$PTwlL7J"
        "CwGi7pQkZ3V#i`#NC)Yvo0fNJm|dH5mx&4*-f09xG-"
        "&#H30@_gUScp29pqap@gAzD#fVPyq8vu6EdwE{bi4^a}So!Y#r_$rSVD3|yeu)B!!wTxWZ>-"
        "hD~C8h>sUa<;`MC^en`p`Gim_~Qoy~0PWX4|M}I*NXPDSgw3O2k`amFo+7XEj9Hv`=6N6(%Ao!A%GxJdnZfo!`uS4&$1DPzNcb`s"
        "E4DayDZlxwip)armoSBKICVc*z^ZMl6@9(g<u+wB${Ud3}5Jfv=QiXC3d&DMFzXqwFh$8rN*V;`P}Q{UrfIP%H+g2DCZ2b_v6jCG"
        "<0X*OySJaxm|iM4X8a>ng^?CoAW_hi-S@gB6=F_Rfc-"
        "2)jJzBFx3%z&Nc}2_<0lG!;oY$xJdm!YMFRl<Z2$l&KiE&$!*PL8~f<BqT?33S>jI*E3_zzL3*~c7fa$W4DFeWOABJPL+9=TXO&S"
        "<kj~rLTr}#t#ibk<1~bpY8%QmN@lBB0mBYrr#Wp|7;GU)gC$ll`IJXkq~)c4TWWAby$Y~{?AyOftaQ3n2xlZ(tOJ~WXAmN7ByE5h"
        "<YXxJc~6G&RhcrmY1jcFIT)?L<Y1d#FpAC4hyeA#Noo7KYq~682n5wM?k)^ZeYVHO{F^SI(QSZ$Mn%2MzywzO$2Bh&n;vL=<Xu(e"
        "oLiK&>~oxC>)M$+H6EUZpw+ygcs~r-$l~-"
        "W!?>cS4}8E^wD~^f`;Tbb0XGe+hWx0d?G!Cy=nY=T`EAptri_~6d;umnQZ_8dr}6%QZ^&l}y6z20H}7$Beiv65%c7>k+$`E1ly7("
        "KP1!Z4WR{c8f-Xxi47^7r3jwxdJ6`;_#w;3WQAZCX4bWhGpe4$A*$8tKiz0$8(u-oCfw{1;Qew;>nFe)P+pXCFo3w~3-"
        "|g<XqvMr>@+jzAwthbq4-emW1hv4*(isf4w})UCJj>A=$;-1AIeb7#zL@yFgOY@S?6#gnOW+Pn2|sT^2rEu2-"
        "0UmM3WA&%t24ZKu>giPBDAX5VqFpEaTWO5dr@l(T#9LDwI@g;S}nbMVs<M}0(xjK8f3U97UF>%$(``kP9|6aa21%3?bW4-"
        "aok|>a=#jpo#|zy^@;e<(?Tf^Ywk+!9hG2xYvIDd9pbnFZS<gW7ULXk_Xa|}MzcXdU!aFUksESMKkR1@zA%c}MS!gG`<Ks@iv|l@"
        "TG+G2UCc5sE<$MP*yb}F5B*4l0r6&)U3m&zukX~#0d##V=pBf*!e6C7?A@Cq_Q3ox{wswd=)qQH)^%PX*eJ6IQy{a!$1ExMX3+6x"
        "!66or0%<)QNbY6>eS~A@c(K?3-A|CtR$bs*Ensqpee*|hRr|A3S@o{i_q;DdX)>h76OpJW-"
        "_}6b?&~TNKqs(KR<Q&=D@CaT&dM!Y@qofHXqhR~l(Q%6y3Ou1iJDVU%CABwGwf4QR=I3NVp}KngLJMJ&IAQ7-"
        "m$m=nVW;c>qcR68Dw>YwIe)=Nm#5kgM<;ho-xdBc*qkYM4+ICK$m(5A#h&Kz2P3_@rJ(TQNSoUxS|@6kIQo0azoWy7zOBW8iU_$f"
        "odtzR>y*0H68{*;JS*$4Atl931IFeq?EMDz6`Gw=Wit9qH@#kq;2#-"
        "O4Cq;qJa11ZSFKBm(tR)mWL#CAFDlt(^(RI<Qrtj$w!*!M>S6hHSnG1E3Nr4q&+e*-C-0B{ZrdT;wi5x5FF5I149ZDf`P+A(-"
        "g)sNVMiHlO9W)wC=H!f7D|Ew9FDRzwV^0*#ZtjY&^v0A44n;syur#?eSZXf$RvQJsu_FJu{AMW(`a&S~8KP)=*L)+JcF7SM0Z#{U"
        "HAoJki&sSkV`8mBsU%i~J(aM`=`WLlw_&>{Oe`UuYhm=)q6>gM(*rzF5<FnmG=EigI3F3GI!2ZK=B?!;{NOT?*Shl`k3ht>a5r!h"
        "h!&RthA?0ozOn$;&!iW`*Ym@^gFz^7xRA-40jK;c7yrIM*G-=l%{-SAjwUzex-iNGQ8?!~?9ZqA3BZ_;S^30=BMYp<PcwtTWI^r0"
        "sMXfJ#oSC*#xp$|qfvWkzX{UNRMk4$@P5Eq~Yb>HKameoyoqmjWB80h9|tZqI{W2cp+;snV#&cbVFwX-"
        "k8}xCGzBb7<@(7RR6y&WoB%mCxMd7%A8JI?EGg9S##91R~O__%4HKF4>d1&e&u$`G$Gt^dyYo$=}1`TNYWk2B!f$hB92l4l!QIGT"
        "3*Al0R12Iv9E!jABD+zfp?#6~^4Uoh}KH>)3dX>ze9o<A$D55MhBSO1r0QEu9I|fqCP{EP2ZQo{g5vfIX65#uNDJYiNc(HFcGlVm"
        "w_z$O@-hN|oxR>P%3$;)<OiLj+-&hh`&r>dPG;lp<aKtYgxE+3{Y$ib&+U{2iG}B|2vF1<|E2oOtf3-4|5<k{GqU-"
        "(U;e6rf|C8GN%o_IY%wptETD)17jrDwSge`CbZ=4Uuucx&-*--"
        "!DCmleOvm8RW%3zZGo8S|W>ywPs(yM7V^+#7H>9!U5YawDYRQB~-^y19mC$CQ@X!A4=74?{JvY8p1ims<XR|G(&MsRP$O&MjR-"
        ">Z3BzFB;c)nD$Aj|h;4zB>6R4BsdJ8trxC6<{5P4aNE@M^jKTnyUx?>K%+vAX5LaiY4=X^&v=A8ynz-kJ41`s#BdVyL9mZfeXd^F"
        "=Xz|WG1=06S!WcrQnyb14_?XafBkjrS02NgyIqB+&AJX0(IED6Nyuy@8s~siB8+4-"
        "yvwAiZrONh!5|Sr!IC+SXa0XNKs)6Rx^h3SpeMlbq@b8DhJNZz5I6mF`jqm+t_MtGX{^i3z`GWw_XzzKhNCqVIIL~pj?HYlejQ0N"
        "~JSD)(CHxh%9wtw5$QyvVl>2cLS_R=z*PBOA1_lC0X2vEx&O0#+dH*R}d3!YoEKj^^vAg8O09JgO`y$eS_<W8u0gyMo-vuWac7ho"
        "yh1=NoaWl%(Vgx$f5~^fw^Zs|vS}Y=AJRq=GFTI8a77{@x4v0CplzeBi&8T0svP1mfp<Al#7FdQg9#Gfn7l`>o%_}*+_H-"
        "mjF=_|k3v{u55==EQ%-gh~0=xrITFa~!jqVTj*TR7%(zL-t_+wt~%N{qNQ^A<DM<O!go6eAmUr*6>nMwmJ1~)cMDv76&$P5`hBJ_"
        "qQ^ch8Pjo5oP45nbDD_?qNRO*{ODR51F+OKT3g0KoNv^Fxpz?~!?(|xiJPNuKW-U^vyv)00S>WA-Ji7uc-"
        "#aG}od)|1231LqX=>U8?-"
        "^Vj*;%w#gv^VV>uox8GZ0uWEG3U=Qcw`U)l;BogfO)6tm70CTtp&u60#;rugboFN3X{MgJ8IZDxUaZV%?bbUU(1ve1WRB|Ii3Mz="
        "QTu?(Z!=auQ`Gd@9H^sedeF&8voorn0fx0Uf2EIs?YbRtnFt?(g|y8$ZsboEnC|S{Uy!h+8&^<ea%MUh&*^O>!T6)@HZa{@dAB`k"
        "#z3mtAWcbb*|W0PG<;OR5}`wYfCfYrZ+M*+rvgA2h~7q52I{!kB8#M9c8|AE)INZT_?Ohq@+SswGV4xY~xv;!N5NmPfxe?J8_5Yi"
        "OYfVgX;m#vgajfMiB=*{7kxc&iGo0npUU%jM{RdciM;Z?Sp8=xs}Hq5@E7YzxI2mbz|(5kaa7BM2+ACjJD*P!iX(ywC#_BziS$Pv"
        "TdGJip~*jo%G11SRxsG;{wrP_dwbY4vU!MhPxXc*pB3K8Ob;p1WWj8SHe}hvV2$JU0yRU1P2wQVGw9mqHVLYi-"
        "7u_0R{pfN!50U#=@c-sNs`+HQ)BfcDC4Pe9)s2UN~e&8LeEXB7p_&*2u$rG{Ag^3WkAbkF!&afQ^-"
        ">Wx#I>n^iL)<5o0i{Axuux3tq`mgqTR<2VQEG5H8Rja__||IpLWcpv8RZEi*<v1YnzAA~)5|Nh16AD+BFdh?n+d-"
        "Lwyi)Zw)n`(=(9x%cJ=wdBsxhN~#^Ma2aktC{8vuy`toR>8hT<o|HqPS4+qV5Anc#+VtI5Pk<10e0ejnT;z;qo58KYSlLsiq7>9-"
        "m^WcD*dXo5JHwUQvQ*;Uxqtu9_?bzqi&Ys+aM)t}_5W2Oh8Dx(1&qJ3FT>cCp(F{tIG&H89C~fxEv8)%h7v1zuY@jx9Th(qsUMa9"
        ">Th2i=4Fs;%5gV|p`fKhHk$R{fi#FgBW(W!~?Os8JQO09ki@E`1U-ud;2|>v-"
        "F$ELm7Pd~v=ykHuiLk&OF!o|1wAy#+lp3Li_b1?c0NpKWO`-"
        "$IYc<8M^8T)9$9yr&1eh1M9nD#&{<+T^i{beC9vv9&k8h3@0gL?vv}n(Y`Z#jt7|s>QruyBgOewSN4ok?xQ#c0ym^#KXs2O?`eDs"
        "UAJe=`9iW8?DL2?P+4~2YtNrxF`l{SNAGH5yiX_yh5y`rnxMmOOOXS&iiW3jbeOn<Q`YZWtM>b&9i0JsN`Fem!hLLh9yBrVUdVhC"
        "24e?>Gu!2@9XSr8AmC;r>6S^5i4pQg(J0HF(*BooMBS2VgLTWCZd}`Q|@Qb;<wj{pogNLyj;$+s1zuPr0e45Sy9iB0FnW9lQt2+g"
        "YNcCbqHXlE~*0woWcLlzn)=A7G-L-qzmUD?8<r6l(PYZWJT7j6F#tv&Y><NL8VQ$Myn&Cg_ySx<v^%|7MxPOBL3m9{x_sIlxW$qs"
        "PGO`#e`SieNrz|psg<5LnR!^!N|^lwe$92y<(k@mm)hB$PCFQAFbQOC-"
        "h}BgV75Au@0<>9l>0)IL|(T(poowzZ$Efjf|fap^o!fN~nlbF>6Qq>BS(>K-"
        "+stuXij=p?ABTY4_VqLRu7IOr&?x;t2c0fBsj)6f7MU2t(Kjyp4188PA$?iBlA?)hGgGWwVHOiz8v6BCgJKVAFBC>3X#?P3Qb9E4"
        "<gb`YMmaO%*dW#=hyeik4Yn8p0?<8973?p>b?+I#HHfnd^ikPi&7CSSm|d6+zLTQC&~%c<g~V97hY)%Zfz~Yob=>xz4kNvyt#mR8"
        "6L_!w;Su%y*&E-"
        "kn}XR{55WbXqH;(WKLHLA1cIKm;gQsD^$UsBSDO9g()Q2Yh0ADnOg$>9oBd`gTUhzlp&67JI8}Oi<I}3^WQj9@$N2ds`~$R{X2M?"
        "adHX7sKz3TLO+5xhnPsue1h0OLAL7;kepX0^QcoF2GeC5UiL1NK3X{FFoB~V<bBV&PXk$$z&u}t-"
        "(^UH=tg+&CuC8HutdDN=vwS9XBq(eM}K88b8wSLVZ(+5Wy_;dMcdVZN8`&S6tVkhPC0o6qISO{LVrYR|MvU_$5jy#8KWjK>0R$t3"
        "#Zl4q|G&eT{^Tq;WtcCOD>63{$!CPidKM!wZDZ;&ML#Gseq%z-"
        "lcYL1H>df%gqG)!IJW0fZqd0q{`=?AZb#4G2?lQJi~vRcE9kRPTU&9L+}X_jvZ&@1}(uQ{-"
        "qYIyFo0>`qvFBfV}!+|#=8A>N_*SM2eqUDfGaS0Jv~#e2W(jyDI_Qls7o+hIHeGhU~3eDI`SIMF!eyph{2cxt(lnWsU!bhAr>Tf@"
        "RP-"
        "#@=!WJt;7KY4Wc#|+1q2M!bb3!zbKJ6>jW<waGPzK`emre_O+l)dGLk~4`8<li#KQMFWapq=&w>O78x#|*2Z{GUHO+<M9Lw!58;i"
        "?!+vESV>g<#_CrD$&5Ln#YTewU%HJblxl8+!<^>9#5=b^QT6cx5k%h+RO5SY8CO46@4=|UA2iKUpLB5r}dp;(7T16c4W9cFzky9+"
        "acjCvDTksq-6ypM%1m;^Zg<dofi=rU0hyia?Jg1^F0EXpYB?h)n!~I#G1Q0(R$yi?X?QiE3j$}W8&6AT(-"
        "*a)=SH}%3C_qCamOb?8Zb3hK(_%)KRMKB8ynH6{Z?5-acHS9c);R8W=8N2lkip@-"
        "2(wB~SlEPJK5;O>Q$X0$lulSKV^&+wS?c4wh?gy$E*334Wi4XAzAqenISO-jz_-"
        "+jHE;w%>QF+xHV)y_4y^G6knUhillzNM8Zn)f)G5AOG+_{tNr#xsdG9bN1@V@jtM`cW+*?!=vwC9J6;Xj^7@=d$Fa8_faBF0?|oK"
        "I*MQ)2BHnO-WmiL$Q(owP6XC{U1SL+U$Ar?k4}gaL3cW-Zi)}A(>r~Kd&JQayo*o2cCUG!nnU*H)&aKyK)Q&Mo3(%Uv+X$P8mI=A"
        "Z@oa-"
        "N!NKd+aPu+RNE<1OO`g?A0q*)Guo7W4~q};U23y9!5oe#FTTWM_NEJ{gPd+YZgH}IZ0C5odlC+hd?ld2>0Xjr)t}(?9Pve*<#c{Q"
        "90>!yOQVdhBU#pNt_9J%+^(^fv+p@;B>ks!c&)olM~@P5LZ{oj=SbeR%TnZ7rN9}g_cTuH#1B?s^&RrK)0v$PZS+|kJcm7?MOmc-"
        "zvrqL-Dx`xw21+#`lHRiPUlv~&gm4a0;EGPhN-o4KBYEGv~}qpf}ztiw!Bc$TMsTm)YxCl+RMS;g<U$JOG!lE>o&9-"
        "x}~NVs!}`Rh+K2ma+I4u=>7V}Z=-KHuCREm<92FMZ>dyp7O_>EOxK+&m%F`|v3%gv9iNme-"
        "H3tX{r|IOdStb8CSM(VGS@9xMKNY7H`c)@I-"
        "%BHENC>c!me7l@4>peF3rv%upLWN36OoezH*4%vp~JF5@J5NP{|$WxBHk=hqw}Zi?$$0{d9D1{Y1}2-"
        "Jg#IIbHh)Zj<<QBPNM)d#RTXa7#^HuF!$sIe)apw8{h-d-"
        "D80{c!xA1myXPx9`7W@1DGV`9eBe8&)slyeChV{C*9pftW~9fqk(;fQA`jfeS7+jOy8B)J=>L1RB;kI!8ZNQm?_*N}pXOzG(n9k8"
        "4rFU_qCz&<m<;so7>sUh8Ip+lrV^$4fIf7x$QAum{?8Y3{4<^`P!_d5>^8?sizNMG!hR17Kh48&>R>RQc~r-kg7)jPBrAJ(Ja-"
        "?rCKQZNWhO4)AaB$DN9)FjEb60!5#IbWyf0Kv-gAMN&egYHeEVb??VG*miBevo-kUMPti^%p8dq`sM>w4|f7`u+PS)@RxZsK5dzB"
        "xA(6djHmtQJ*6Zi_FeYQ*)J2LCs~?`sXLU|cZcomWhwzhzE0fR;etJ2p1G!ki0S)Veo*EZD@E11nDge6^TI3`OX<;CiW*8F_MmeB"
        "QF5pCylk3s+23K`0@;t6)CWDtY{pT+FG-VM;W5n3*a2rvbu&>d!Ln{CoutFF4_%zMQ0pkgn)7p@SS)l@Tx|jhDd{*Rr&20VD)e&@"
        "3=^P^%~E(ns4VA1-"
        "+#TWsIz^}_l}dHRu65K?V}=szXbxmW6Th*iNb5!Elz~3z`pPNaho0VeV;ZDEIB|dCi2y0vPC3Y+7p|rIrV8#e%ZIc<sj`vj*m~f*"
        "A)HX-"
        "~Q+S`QQId86a4|Xo3hTSW{<3PQ!6}!~XDJ|2w<ZBa3IS>b&cg9=ox1c+2;d*0^((h}fVye*5Ct`*%OQV#lxEym|lK(d(BziG^i(0"
        "VXnnAu9A3>t$Ku-TG(VRiM)0c;Qy})<f#YWIv|+?J|E?y%!pa4H3J%U5~2l+3!Tw93S29aKSO5#S>=`ZHsg=ULE1KCwyRiX`gbzA"
        "KNMO3vP%M06`Ky!~4YOw4J(sP|zj9dsmxkz4X)vJb>2B#?(Cm=FT2J4%n*9H<9=ys7Z(5n!<F4i3B6l6CfTUTTj0LO1J&*ym3?Y?"
        "q!F3Y2UQjC`KHb0hl)i`o<So^XcPe|7DruwU?7Ad&QG1E&}!lJ~}euJ1oXYx7aIWoO$$%GIP*>Nn&4Yrir7vy@e-"
        "5xrES<iRb9xxETFZX>n-KzV*_&8QQSZrsS?LT@YDKF-BMIt;{v)%UNL@MWMf96QSwLN<aIxbNfnn2-VIEvMR-IRUk>Jl;ez7_Grj"
        "qye&fP!x}>W-AoO-"
        "b3m(y8Z*Z|pXMBk^aOy=Bu34KxR6)RPF=PzYUgXj5_Ep>R^DML$FWyY+d@g8lYXG`XLDU}TR=Q&f%wR)g6FO#q3+;C$&bBb+PXxH"
        "D<HMi$*v1@I_&rX;+)lrj=VY=0na)zP&pI2RQ*xtWWnhwSCQ7f0*78pwj<2~X@_^J9c?zyr60e?ZJgdFn3U+hY-"
        "1N>U;n;k3T{X=N6O&doS;{-"
        "!5=_dh(K}nSX`xsyFb|W+i0FOTPi5jCXq@&gI0aIkFaGzTU1B@MV2L2u7OkCH2Y<m;+2EvWzv6t#p8S3U(r&_a!$Qk**Wb`i^!fg"
        "lh@`>*>mw&hP704=Xga9WtfL$52%38Ir*AEgN4~U`*t7N+)snr6zr79ztBGx${V?lB8dE*_@!%cegKrX)VncQ>`eo6lqC~0^^@3!"
        "`*e^)g|7O?NX7Q&y!KiUtZnJ<OXBUS!XNY&<uMR+truJo)(!#{;HUuHDA9;1MLpbs{wQv~k)2j`mX}mF$oI!|-RAKUvxMTv5DW-"
        "|<4^W6yb@8(futMnQyj59d0pao<n<x&P>s;ahPo0U;&IGN*zosf;)9K^8W~kmo~=Ck9<Y5!(9B(<qWQy_D~R^HgVsAl?zO@n=YN8"
        "~{dn5G+0hk0b_Q}qg$wqTXO0wGqNwKJO&()yHhuNg&C;#RuH(~DOH0p8zL=W%i3&pad}r};8KZnSvWgq0AHm4x<t1LzEN-"
        "t~lohCax@%kfuHW31WTnY6GISZ<+G<^XDhqp>W>rrQ4lXY*!+4!!;6csneO~1Vpa>Iw04S8}YQgJ+S4Zy-"
        "4qv@Jc=ibXJ03@{pq?2V!}CV`+<QyS-"
        "C4!mSH1nB+H{$c+ZI>2dLkPMsltMXy*<w%^uc80W5nCh2&jfup6Iy8G{6IUJm1{#s@hr#_vcvH?CvW&!rKnx{Fep(%L4y@d4WHgj"
        "7B5&;Fopou5&@MC>2F~;#)w)#|0JE$<&Vcb!%9w_zTLzz}o)C%uIcrrC)f_t9z?;0!GB!moj$&3q*<wrF8_BltF=J`p0o^NwVn)Z"
        "(D=yd|O(Z<#()Z^RRAjF4%XUtQCV*yq#?(x*5F09~$oD@$4U?l!+ek$Ra#%7NIW=J>ea-cz)3?mJaQ@XdM8VS_O&Dn20!j4iPSWn"
        "mo5sUdYToogrQzajiL{E}?VHcP}<loinT{FlixTZj{Ws2Zpge2l*<_(V#|`=}2Adj8Qq~x)excg)~WLGd3*Bf)A%{l5_ll3{kRY^"
        "SX13$ks|^@#`(RyM6GlcNg7aMNiIV3sSJW+bpFFiwkUWPhhRnHMU+dci1D@4YRSzS{8e#WjVv@YtI@BfX;;KaI`j61B=#Fzbl#~v"
        "++S%TRQQ0r;PRRnS`=0VM{~NTPyX5b;UD#|JFiFw511qagQu#Hu+AWKZVq)r?h#a#w!AoiKg0gJVYjgI{AH<;#6!9N>9&Ao@}7(l"
        "YaVnDVF}6chS;{&st-"
        "$zTZRH?Kix2rgGPmaq(qS#%(bt`pTjZ!cY4ZYJ71@*=iz=HjJIbgXsH^QW%eYv7eey?S_@dh%d@4ak3a$YJ|kDt2M*4_eIpJ&c4W"
        "_N3DGki>^}a`bOB-!B0yA_pB}pExlq|z2tIHPya7AQu|&"
    ),
    "hscan": (
        "c-rk+ZExE+68>JG|A9z8SV|SwcA9S1fDYG}O}D^qc9AS@KNyBUOO(y6-cFJeJ4OHd&5#sDQF4-"
        "YcY!<LilUXJ;k<G9%y7t~qoW^BZ$g@pqRc~9RXkghc$<Ya&$DWLbaeFm+4E;SE%LG^*)}cqgjOUg^t+VSMUvMEUyX}BKH>KwsfP&"
        "GC6B6tg-Fso57C8rl9v@yrScJVh|4@BRmh9|IF4(gM!{pX;c@K|9ykAr7dS)>3iFaR^c7!cGzkRlh9@gl!b6qjdA$ixox0qn9_Uh"
        "i3No6(CwjZifu1KdD?J+hv90Q$%y|&80)~_{TQiRo`5i0u_#ojKye#8j7Oc`dVp5mq&mtC+HQj=fX%>WerPw;LM`Dl)x~0@3C>1U"
        "6FKJRScmDhtfj=-"
        "4?@C5%Mm{JGU>D35Ba#UTSqA^Hh*0oKy2zGeIZ*!lrd(Hw6aI;blNIAEqGd_<?wnkf?B83Kg)qDL#vn{B%6CExiPZqks+L$7K%6+"
        "qX$_<|HN|leiP$Zt7*H)~MjW0Y(Z?)eJCesrf6ds?K|WoxdRu0u5v^YH15VMdlMm~V_&AE~3g;B`B&X6-"
        "?#aoFyd=(yoJbqsx9gG@9oiQU(MtY2mNRGQC`9U@bk?CFErm7k*@z65u@WXY5HX+9!ijk5PmKWhpeTr1vJIh{=^uHup(%L-"
        "jv%wi?3_3o@n%fn%|B_fO*1~;)_Iv{bv_QM`#^wz=yoFC^E4H>KtmY)_{8*<^7toS)hv3a^qt!VSo5UTznH6gR_2xC%+6nX<P0W-"
        "920m}B=Cxa8=7S-sfcs!kq;hu+cE2yT=KH2$%ZEJ98`@+os+i-"
        "4L`~6ab6Nmct*rZcBZon*Dyp%S@0$3vY0}+fLl&DOh3?hWCDMeDB?W}VQE$j`FDliurSXeC0&S>^0(}ecn{4xsKZ1^Rj+2}=cjzR"
        "m@I`7A2X$d<lYIhbZ{D+#)%nCnSZpYU`mc$8IZ?B2OBJ6;Iqb|C8;;em^Cp($15n<k70Kpz{M~pE<tQi#%NC#^R-"
        "LnU`CAk2zWLi5%Ji+$`MuLq*hONJPQ(*iC@Mdc~cb0p7fI%`2G5l2%!%mqteQPy=EwhXiXdMiZ!1p^$i|^Ph|ihO^ZY<pgy<>Aj5"
        ">F@bGGzLQ+CNyuADk5(~VoEY6QZvWP1Qf!5$hA^Px8t$ZN!B^u<c@;m`P(&?IlXbOu(qV4>pKxs5-"
        "(pSANSeIDMZ*Ks+kMcVLpVsx}>#~H1e3&v!IYy%E0#6{n8|c%-#2C=DfXItbSsy-"
        "cLA8qx-k|p|{kLnnJ9sf&nu{B11VQcVxLeFU!C+}yCb-LQ0fdWhlroND)2yU1E+X}avE`ZpEB(!)pq1UiEI74c=4uvKZim6Knr%k"
        "IO`CVsG$L4pd0K31bipkI5yL<>#Djyf);n5+wS3Gn32|k(YHG|)^<z^*xmMscMhnk1K0N90JupJK1Lq}MN@oGiCJrWb$KdNO=BD?"
        "S!%k<}t>JaLsvK2?I<`QM&>gS*iR)UqK*pzz6UC_nyLV1cj^Z+Pp$NdJb*!i@O8DUn>;Z5R69(C0eAF#j3>}7neP$qg45H5&r0zK"
        "MXJhfSx(ndPWc*s8J3OA@`zv4M_$zR$|7!BR&VdzzC#$5jm9bLff>N0lhru_f%1#Jc1ZIYeC-{GouV;?p;x^tqj)B^NL5Nw-"
        "l}I>|KYbZj4F$lg@GOX}3K}xUnVt&W-IF-"
        "<h%<d*{JPR?2S%yTROC&}N*)Kl>4otiGC5uExC7X92$$r5r(%)0{NJhdjkN$KDacj$2q_YYzXGigrljl54bPLoO51*At`Q4h`p6O"
        "<#*@Rlur80|zk4U~C~k#j()o7G0?<-kT>jU?)k7!*1?w@Etnd{}*0c9e03sHV>mP1@{doNYIX3<bfmWah9jxzU!C3`ok*-"
        "bWR;6HU#e|$3Qe8s#m8zDb(-NyGQi!@tmYl=;ED+Z(HKh^)--"
        "xxuJ;|?f&Fxn833RpD{PIznflUl3TqQSJsS388Beah=@06YD#kQ2yV7~f;g>|ERzq*pn>j}5K@|Z0i5tD3S5e61&08s(Ew=VMm_i"
        "6!>q32rDH7plm5uITs5^c|6kTtX(ECemxO-7A&%}7~>n;0#jmq$CRfb1Z(-"
        "lItUac#}LcIdJj8bH;A^8lp@>^GlyflUlCidTfCMQv<Vj4!z7A_%mQX7g=U2U(t>jSdpH#P`X>f?Zi?#7Do=GK2bJE}2V)*s&6ut"
        "It}qFGm&<LRi*3Yj?I%<F_C@!g6JmizX5Teft^XZl~G^cN@S4N4mo&Q-k?_>XH2nsP@%3<(UJMgHM1?h}+Yaa0i4tB;0Aj-5{Y#-"
        "PhMxVo3I$8tvLtJ{j&Qs=~OZr`M{X#h4DkRUejl(GQkD=%Z8nema-j=DV2{Fp?6VY${Vyv&6Hm(-"
        "oRto6Ugkejh$Iu<C0oP0qha)F3O>E{V5EA{DpAD10&27_7NZ%qTzyS7?tx*G^#9STVHkks%8(U)6cZ`cT+gWSXCn)|oE`=>oyR-"
        "~*`Pkw6-SKk_xOa+;pB)9B_kljV-%gw{0^A-"
        "RIWQ&!R<j(lf=9i?+KEo*?1BC2hQjly~{pUy6pE|et3OBYvt=l4IzgfJ+H$YeY@1lWs?@%ozqu`0j}gD}}PIMv;8i0_^yq6wklJ&"
        "0IN*n>1e0hu?WJi=4DW<kVD)lEX;UdR3ax&D7J*NY<;3_E@PV&cI*qbUz1$oj8EGZotubx`NQMX=(v|693jHzALqZi~3O1N3CPeG"
        "1kN;_RQo*vIkpPhsoV;c5>6|5ot(D}nDpEK}x7CAH^)A)-"
        "ugMLpO*fA#X*?QFg`Xjv16Vs8dhDBAL}4J=+YTXHiq8$&vw2~O0#N%Ha~nVI5#Em=WJy}1`P>MgIfG~xHwksg?5hzKb7fUdFB>lS"
        "jIoqWo=L5sOlH$|kTP5UHE=AHL1^!ugWHMWv)q7DkS1R87)bULuhfnYDrc}a85dpXqV_?o_K(_&DFqf2FRnjEdLYtlWLb_Ncb9VO"
        "V+Cn7CMsTT%s!PU_cPQ@fZ3ER9e>xAPc{^@v3o>7ae!a$*VLq3!|y1`R_jGY<o+_z~4j+l`$zpKoq!y0rel6ovhJMdq{?_m{B0YW"
        "ZJ&DNqPtZopD8vzcG!QTzz>rE{edgU(wxTj02!u)%*eDxP=`4fw4rHDQ;6=h3HM0ZoAO|@N#dS_i>NKcmHr1YksSp|C|un(=kJ!R"
        "4|?XVg1wA$o%j!2jgrGP#D2(WiQ@(9}jzSjOiMb&P%m@IW5!>`k6`}5xVHCf{F6cOu>=(f&BhKG2H>KoEksVAg&!O>#kkz3e<`4j"
        "w~?f~=8|FoKJ*B=#wnt;~{b4JSd&h@oVKF=rPS^UKvjbqdc)_9Xt+yMr=gsxcPjDF`4%z=}p?r4B&uc+SR{;x-"
        "c3JNc*Xz75*5bb=51xa3VFj&B_Kia^rd)Oq?WF%@Waf8Q)T=fdlVL>B&ys#qdHz-?5j6+e0{)FhL;>c8T)(^3X3Ri<@EWT-nip;<"
        "C8q~FD4d*JsqXCLt!c&u@Ibr+RX-@;w<M-N+h9?8M3m{zZpkRE?PG6G~jm%-Utx(r_+PLOdf%-"
        "QJbcZ~2*7wk~=b<y{p%{ctK{nJ+E65tCtv~9D7|CpMhe6}d?8CPx!V4DR1_~#2PK)b0dP`uvB+R}zfVT#^E6JHas<2*7XWv2jgQZ"
        "gM6TDQPmfa>h_hg@)W!qE~Esh9_*>|GmTH1YgrL*d`{{?4yn4U%Hjs;L-"
        "{+!GvCrQ4ZPE>Q*6^T1TTcBhS#nJ?+Aq*bUPp0451SU%!Apt|nf)_ZA({s0(lnNvy<d!B|A%))NJR&d>{QbD3{yUwTd!K~FE8C*3"
        "jDc|L@$szVgXWKn$$0d0C(6j@?8WrUNZmg1^|chWzpM8N^GBzro!Fg@?dZMaFdH}{DJ?(YB=8iLl&q@;ZbuY(#cO^G^$mOzhksax"
        "o7>@CAWIi)nVLrn(Q45hyi~_9vJ=?-<^ul1-"
        "yy#E9px@MqL#|5RuumH<GpFFo`fztT=^C|kBmx8YUmYM&{$8rw!njCS&oMrj&<fU`rPF{$Ml4njO_W{?fu6!ar>icNUZiq@a1HDJ"
        "}@5Q_PoafcJujN0h{E?k9Y_{Nkb;XXawmK6cq$b4E<x|dX|xp*Mw81&Fv>Su%)t`W#ktNy`AVSWChYEZxTm)719AGTKW2{%}bz3l"
        "TQXn_DFJy{e`LZ6$?TKEB?;HnRnLBWPD+bN*HVnBPYp2f3z*I7b21qruxP=7plbqS>Hah6ezkaq(>lPbn|PBSHE!Qrr*ddZgU~v4"
        "eVt06|auIn3ATi<&Qg2!}_CDp5P?zLGSOsKkuw~%!;+BbYVZ7Z)?ZxeiKz+x80`;DzEAMwVl4>c25Rd{hC6#ld*hTbuK?SYzGeFD"
        "b_}bE8BzWPD;=j6&Rz$!Cx?o<f(MZ>P;Sr0yx4NYh+%-eN+Llg?KM-"
        "7(H~?Oo_BwULP}HugB{?8ri9xrP&UA_?8)FXGA28YQieuSeHQG8iyx2BYELtH~P98RUZ#I!D`Ck*^!hED2941TJc^hd}-"
        "cCpd2FcUzksS0C*^!qn{;*bX(PA#SR&8)VYwA1-cIjkJ5qa+zecA9{rW9+9q|c-=!#B!vx<qWs3y^W3Y|Z8R7{~(X4kt3=ph6-"
        "qTu{1(ST$>M7e4iIJfwuFRT;fhKKE=H&}_=>FIN6YQh94PGM^e%NMrCDeeJR<ku=hb7g%uF#C^+?sd9TTk8F3*@~$|LhE}%_sZ3j"
        "9BT<EJkuJ8mpx)jk*m#aapIiKmp$iZZR>tL0!w6P#H2DKroKgso|i}i|KqRB?2z}(fc~beaWyn)7BoxlyyQ5SYqKKGHI{VHjp-"
        "p4xrk;a}i6oht@vS@=t}cy$Dz=Y*j}D0~fGK9bA0FP1<&_32pK=(zgOzV_9=o^mp@$sIXw*Lh4iCO`~ugAI5tN>JA{z(6zM~^qzI"
        "^8VrqgyVtm$3?7|(T;B*D-v$O>0R+w59fAYtBW1z?B!IOuu!1_PN3!`3%@D!n"
    ),
    "kidney_aec_202603": (
        "c-rkeYjfi^lHc_!aLiVjoX{Le6Hm>lj=6F?_RQR7JTBYG?#1!3XbG}-rbvdQ{1~n7w_kSyAOSvPd6Roo-KmnX1avnVjYgv%AP3Ph"
        "OY@TDd?3Gt`E3^FMfH0bmW#nWPnT>Kmb{FXoXH;bTl5ywD0}jvq~v*)#?Ym}MAK+C&?shcTE@|hmp$P(78WdvOF^gZRRnJHXjWu="
        "s_@f7eM?r$><O`xOf_X;GJ^)BI1>ajn1Dv@teTJHI8{Ay8qQ#JmlfeXhc7x*5#1(XtRNRre8Y2qyjZ}`OOV0I$)(F~c^N!J$t-"
        ";s2;6C!%%fX{bws~ZPhJVg9x(VIay}1${b0$<5QH1LRmU<c?x-"
        "tIBhd<=<KSkQ&bT{p1kFV;4U;C(gy*mp0rj~ZZAdXdfj2s*X`JTEFazR)!TaNr!=HY<41W0Z<?P~LgTMZG@mJ={AyA%n5X_^P2Z7"
        "_|yh!7FZaW@m5l>3_hd^I5Yg*h}>ic*&emguEz7>r8{3eK&;VlntqBJX>f^l~6`Qz7Mm0hN1Zx0TFuq=793d<-"
        ">yvMkBw16|k*ZJXJg7=@!fv4TL!O>B`7JfMo9|HUlOk*A<d}h_tJ-"
        "j^n7<~EZ^6aNeN&oMeGmWHve7r!RUw$~s`7*udv)~Zt72HL$gg*tt>>x-"
        "|?dLmvnunK{$EQCXUVi#=8eATp9iImmpFVs(3_kt>oDIEq{B47E9pU2W`1J7n)0c~C7Hj$!KUOreCd>jQUxqj+e`~pQXB^(}82dk"
        "jtD=ZN#|Ti2$Uc76VFiiBR{(#WYYlHz@SqGco&(}#zT%>5<90Keq!PG`%hQxqjJBXO{(wpugjp6x)9@zdRu4BPVO-Sgz;5^kKTMk"
        "yJNj~detbk<ztKK}8J{y8;7M0r*jahbr}*`aSciEYE;;OdF^I4k;`=?%Z_<MM8V@=pVlhxt@^=+D&aOqq7gBd71Y%i2VGogeV37Z"
        "|;>q-Bf5tPO%y0*31Xw*IR;CQZetse*1w%#2ExWipx%3RW-@q0^05};6rkP4H#W-"
        "7pVkLaEelS#SOR+GEdMy{U@4zAO1h_BwZ4bEY1%&azTY<1hGqutk2+#|_aA?_z@GBT=eBcP=<$PZDBV#XMkA~y#1o*sQTgI9f2<#"
        "07?38VTEm8tL_-"
        "?4+=hLmUr#NTtfXi9mJOm1%1DbkCssK*#;vfqWTt>jK7GaWz#ao7t!EG2N!R!WmkRfU>AJX6<d{SeeLJ3mQe2dz$s9>jQ!Y9@Cn+"
        "I`<Mry?N1=kKN!F;OqFCx+ZZD;?Mu>gy(_h}FZ4DuUbUaBORWOQ3lTUmPZ@1W6|Ae7kT0dH&Kag_?X43_^UN($OcLhMC3F9V<k*t"
        ")ZT_`aawqC)!cKcE-"
        ">yl_d!eEbJok!c8LOqu!?oQ+9Tu4deF6vHQKU$XEKqzL{&ZTTK}bh{7}nCp6_`0?hVFMvZ<_02&ZBRE9ln%R{^J1s^nZy0y``hum"
        "UH8xO^%!`V71e4-"
        "+Kf1H)FyAwwg^s33AfG@80%ix}@yE>G2g3pSV%u}nf+_02meJX(S8ObJ;9!Hd(1z(1K^)z2doCA9*>z>{`Ak55<$0<=C5W_S&<Iq"
        "z_9u0~{!|EFu>eM|M%M!Kbc@{xms#9}JG_>BUOYa!_z&^vvd1#PDdBp#+P`v6Y0i+BAeudbkz{$8+;Uqou8AIIjZLUja@4DdTVY;"
        "TaKAi0>Pm@UY>8X1JO&v>Lw!s8v<84<sayopCGKwl{0gk;YcPL9v<;)-Q0?-{J_NHL<hrfk5VC0>xD_G9qJySFg8d-"
        "vTRlmG*ef50t4FD<B+Lay8UYKaOtx-4fy2t0llLv*d6+Z7-"
        "4gc4a<#N8^lChb>ly<&_4pr1_LY*N_<>%+AJA4@!~C?t2KsuC8UQjb%`5gnUsCE3jDf9;?P@h5TDV^-9l9I(L-"
        "=>(pE?>#?G0fl1rU~f16}Nuy7@l)#{*OP7o<KgZnIUJgrKPROAW7;i9do+L8DwSq_B$%{+qv7&i7n4<rusMnoSqM7cW<N(lH;mYc"
        "oLrg`QwggmWH{n`+;OF%X;R#0uE7q8VhH#o?6OmIu4m$}FwrTWX0{h2HoY8v{HS#6Z(1n3h3Q1VAsE8H0g5>1X@Ym28l!sX!xe--"
        "<F!U=YFr8v;Rc>5whP!~iUO5+=q7owQ1jFzw!Ij_@_gY)6*>0sRBG8Wz#VpI)3k!0vRY8Ah_?6fsm6vs2;XG%R>wKcxBHO`6^X81"
        "zeq1OSYl%}{G+0e_sr5}2TpkKa~NPG<Vbf*3DCNB7fvl<v<UijrMEWxW1%x?0{~biu4o7`SS`r_t+aX$dR0j0$psfhNS%WI*vR+D"
        "<hbNE(~)l!Ru%uLmuwWl~(LwUBf)A%{`H+0Q7><2+BnhFPEGSQTYS1j(+A32!EBEh%o4LqYNcJ_Q3;@{gK~w+u78p<2#0u$eu{Wu"
        "zuUaJZ#ATD`C+L6SB&*#@Q|Rt?u|@5440LBlrN`>^FKYcw@&2@m>8O53Hlgn?NAa1X(y<B#?%2bWG73{qI|b7rKYGRZNEay!SwLj"
        "*=56S4s#=dyd8Yk)(l?+I}ZAqgRhUsAPeV70lb2FF`v8Ezd7P~ZY6MmKts0IJ!LAa#QbP>}x&PoNx+w4B50FS97-"
        "wsuHZ|8_QfwXRV&!5RQLII<u`jlWF*U^){U`=2Ml-~QxifjA(3`f&OMg6qSJW6>-)5whrHrHcbnfnbeuGNicq6YN%F5V;j-"
        "rm!~XixuFt^I1Jb+&;98=?m;E$_B~BbUbipb6j0VKf^L%D{0w<_0(kvTdSNdYy!f7!qH|*_Cj)lgH@c>hM>}@E^a$HRzy0O?Pb^y"
        "c{j&3*XsPJI5UrJxo_%ja7OhFCHze-"
        "f0FZCaL&WHCy~;XN0BKs%$&9pE3%VZJRbrA(W+(MRL{Hr&a6wCN4F8pk}@FAB|o%Tc?zRBh3338PbCol0n!dRs<wQaB7H_o{ax<%"
        "Fw^gBRmPx82Ce!3ADX+VLNZ-"
        "b+zFpOEW$hwpKNVwIO==@&meXxUjuqXl!%U(88>cqK>kdW#gGDjUDMjIS;|E+fXqjoCS03Cphu!v1K|8&!E<iYUhY`UF8x*quPAa"
        "9tK~AxpKMVEA^y}6itMA|q^>JXN^S96n;7EDm7d2CP0qF9YZqP-"
        "edQRw3|lVp+OUh3oCU&L6E<|sUhwFG92EijFqTsloeChl=Q*0xDq-6PSd_D^K$}4)j}89Fji_hnXc=^!Rf-vG3}(5C6$<uGY&3r7"
        "=qO#xccB_c&52Y9(9wY3bZJC@$!FLJB`kNTM-"
        "A{3mRObq>rRnCrb3&KFdmp_+oA@kwA$v<L9}M8+Gnz2*9yh@J$s3f;Jy+pZEAA*&@G}SU%Ju~O=HJ(ZW`TcNKs~aI^{(Xz(BH^h@"
        "OS8ntA#lE0xmC6-V5v%vNRn`=tkQ0TFqaM=7{<en-"
        "8k>ZUN;f2LyfUwo?5d_k!na*(pB;Yjt$DkmG{WJX6$w`)&uS{Qo$tcR~*sOU%cp4+fyq$GWQ4ES^wWL#87@cOz^wZPXh$8+m!oly"
        "`_PpJAK{r<{5*Vx)v8VRXe<?U1zO+^;$Zms83bcJzD{KSN8su#o}YfZ5>=-pbH4=yK1-"
        "m@$**LsWzWNr3pFKRd*e~gN<kbX2Cf0b6$+f!~fbr7}flSa^a<`mJiX@Z*8F#ehnP4sxP3rZ@(#EeBvNJDKVhITa|h3-"
        "s*QZ@Q6<N^@W+G_NG{V2=5=8kUwAJ#kb5;x`rZtss+{N2SRBN9Iqyh#PnxA<ihCRO^`=fg|n5_XTifUV;=ELR1U;O`>X)R?0}m9Z"
        "Ld98UuTYjnJNCDlbg!3Rcch}qz+Ifxe|9|q0vt~Z<b|BCs)Trr!Pq^aUOjAB0PRY4(KSbT7|Q!rd;*XHRenR!-"
        "A3*|j1UtWQdM;2Xi!jwOiY|3z%6)poK@7C}gm!j-"
        "S>l9qe?`jaxFHd!0Ysmwc#O)Lgg&0tAFtF{=&;&3}IOUu+2)I;}^~Xj&L<VaDRQS8#)Kw5zVQeyj6_SZ-U|#F#A%fNmj4w!-"
        "G^XWt`HKAkui>ESKWdPcAJ=U7pA4UMPu{!}XE&s3Dw;T-Q8pxR`l-T1Q66c5UG)Ut6cb(xcO2rc!4^lr*6B5!@rIM4uAwdO(8P-"
        "$qJHHvI#`1APQRict~=H`p^>U3+%V-L<e3=LbUQ$T6}Xr0g83?r6|~zK^HO9qc`i)>-"
        "9*r78StRjfY+Je&Qqs8Xttd?`Z2Tp)X|Tb?Wc~u69BS%V-RF>-"
        "=Ls+w~~ma+nYunNLx>_Hd?&PdDg#0Uvi{%Zr2?TH8htpIxwyfY8i>G7O{o2G16@9R7k%s^SmFkn-"
        "8;_|I+8X?BKWT?6dq%f93D<Rd#fsHID^LKrt^3r=`<pZurY8GnB8A;0803@)*#{#{phB_B-"
        "xD8#xZYaUr>x^YHs_7XQ9vNa62r?N^}0(K0H1y4ABaKSAYmnE>l5v};nS!U3YKU-zoR0Y#4fdJ7)TumEclhU2Irazk1jv33Ju*As"
        "3TCea*ROclawW4u=9l$KrSG!=@5S6yL82YHqA;`dNhG-|5Z%zF=>=LxzXcD2HeSEd2L*ybQQsS1{8;T6^ye-"
        "z6qMVsX)^ttsTEGqae>zF4tt%$Rc#ecFTa<%AM?pj-Dm1$IDz7Q*Q2ryymn!6_Op>bcB&*#w;&1a#GJ%6q1_=b}D`>LzC1x(@w?M"
        "ez0t1h6bMogU60Ch!|$QZ}P++u8>Rav2{y=b;~jh80s7P67h$gBpYO%~ACK$S@ZXjj}x*L$nm=Fsg}S|wHGQE}&739r9QOVA)NRJ"
        "1v(^<%H-N3AXb73>*l0X`mJR+FjOWp?1=01gTUM7Sg;SkR=}Fb-"
        "}zE|}jiFa42r#fN5Z1VSp^r_zag^$EFm4H~wl#`<jR65Y==wthqszBmnQ0R41C0IQ=0y~E<0g8J5HBeQni)(yex)w;r(48NOgUa#"
        "x5p@!?m>;xLhmaXeGkht|yQUmNwl<%gS1a4Aah?BmfB-PCPrp}ra%tC*eH*PZuSEpFTnDX70-"
        "ESeOY=5gmo{UBnp!v~G$Yz0{`Pm`r)ITqXK=qGaW~k(e51(8=125_hHHz?_d$6{|yfuZ&k?tYV9r@Py^|LHJ1Zo$c`ql<C^#*UyW"
        "hFbNKWKHs0^25dDLqMl3n%P{ACHGayug3QdB(~GXXIEhI0~cW*6UZgbtjLex7}Sf>7JeGjBhWbY^6}U;XS359_*$Q-"
        "&1Vq!B)<FZ~3JMuhaS8rop|jz%x3$B{bMwhkF8q7bx-"
        "8q~@;4%~s5SQU|^65b|cwuzrOXX83HBEi~EOFNRNtz?C~^c&1%Gn*tE`R~u9eH<YT)0{E)>qNaT9@!xDUYxCMBc<z+DRh}qcW~Fv"
        "4RPB;u!|Uz$GOb^eA85y`pv*ITwhjiVS`QKV3cNzPK$(P+E)I?Ap4{4w=~aG~WbihtE=jE@zHM%b%!%v=GRE{87g5v@TDXX_0~Ia"
        "O#8#10XRceGpYowvJ>#g_?C)NyJ6$*E9moT-"
        "{Y>jio@A}4iSY$JoxDbdcF_GBtG)qQ1R(fQDBo2jG}}5G4Tr(m<(Jd*U^oi?WfceGxBKHa`{&0$jNXplyxoU>yar8hqPRP%eha*j"
        "YL}edzoDY5w$o^s(KEv+vZKZ@qn8acerA}lBs*>lGk#_m8s&Omy593Np4BBa^-"
        "JYZWer_k(yUgX<F;pGC1D9Q^6*QiUtvrB#wuoQD97st9~E9}%l!Il;_vpcoEGw#@GS8<yYz(xn;T$Yt@FLnko7_eU6ogYsOsyJ4$"
        "pnwU#Oa>v9xeiD|VDoq8?xrkK2trhLJNvBV0dv^oJI))P21dR<K$k`i6}?aL>v-"
        "gmC`+gwdHn*AK)Ky@$s~t*edlfg`MCiAl)<R)6726s!>VYT{1Zc8#jRlNkH*!x1~joJbgVRiWCIw>CgW&249PHJq1lt>+ky>gS#9"
        "Dp}_=3T4%mOz7M41k~<#QH{xNS@|~}`ydY&_HiGLdSmxb!+5&FH+0PUDqkjp&9Pc9upsqJwvVG=>}H?5joCGm?>rDSp0AueL}0GO"
        "2{i`onR7U`rd&4(hGpw93Ff6f1H;_@->{PiA3=h+MU)8d_!G0E+a%5D7G5-)acV3~JV6qH;l#Wnn}3J49IpW)v`y-"
        "{M__x{aaY|5y7MT>Qs7MOy&%|guGOvhrtuc{3a_fcwm=y|zaKY+Z3vvqSu~f1AiODBYX$L~qU1CL`xw-yFd{VA*w6subk68-"
        "WNuo3ore5#_{iP^s91Y)`NXIb@9n&4z$+vepn;{B=t65~0mLQjKj^Mx64J{J(gbY5R4S#&nTlJc2BJuq%ZcKjtBQ!2>~-"
        "^biyA{*7XF?M)kD%@98cr4pv%P)#$`jWcV!y;BfdX&+BqEQgYWNu&yK#Fo_sq0eEfdGj_QRJX%=x08$t)QvvRew)Xn9Fc+6*@J4X"
        "A<;rZ#O(+?9SlEVd?rP_sW2;k;mRN|*5-44+oc<c<{4Rlm_*k-"
        "N}pNGJwH8C)XORren;9G18d*re=E_S0FeS37^vbPQyF$f{T2;AZz4xeEDT796C1tbnT(2JJEB7NwvTT}4bebccAmheqQrN0-"
        "h%5=|V_sb~3M`7<lSN-w9fzvqzJ~&oJWz6k8ZH<KSX7*`{=QTcst)HFF-"
        "gCT<8ocB7|ATlQ89cY_GPy|G_=}cmcM!R8#i`J4;&N}1mU|86I#?MStQ_7HX}p4yuio+vBs`jF>NjC-OFmGOd`Zn^#RByC4$ExYj"
        "1uWwCWf6`N8rT1z?ZEk;GBC#Y40^<(`q@3tOtl5gG0l!H4Zm=#FH*;<D0#LFu$#iMry%WLcRCPJIrraONa^2sK?eGvM;82MEO47("
        "mV`0AL*nM4}Xa25gphmL5DHvZ7%ed%$ll&^}I&L({JAj4!v+TL&k}<nkg_klvva<IpI{WSb(?BgAz8m*16zu=3D2hgyAJGs;v;Gp"
        "-})6k4TfBL#a#g1hIO6R$XxGXAij|oPByz6~0$s9<Jij$AnTx5!u;F3Ie=+aQA!2$Wrrm62jbXq0%f!Zzm^KdS5+_hL;<By82WTu"
        "RP)&5_Ng+88Yv7(*jMA_9;`iPf5r9h5mXfZ3IoN#YiuCT@`7OSqP~gHK^Z~M106c55U*^_?5N^Y83T=i5E$6mFM&t<sM9nvi*IGl"
        "Vh%CrK(xrjO3uk<dBq$sHmxEg)sL8`>1o6m+vdd5U*ZEqr^~`;fj(&7@f{x4;D0w1f6R2Qm^N?R|nTg`o}S@7||rU!-"
        "l833`&c^PS-"
        "k&S0T3VPE~gS^bTf?9Hb<n8}h+m03KvO37a6mJB3ydEb;9gfhCMQ!l9p}BdHfN!7q>0aH2<3)GTD%xIQM{!=rU8aa9nF+G_J*rk3"
        "Mp9W7GuZ^Z(DtWY!hv>x_PyCk|r4ZObcXJb;$!jpK>++0>JjHgFv?pt0&MDI{?HpV;jj5pNHn~x=}cOp5t-"
        "85gf*L<bjGjuwx+qK`hXI3{_iwlz4;j$4ph`VW%sVh}^o?hf*0J@Im7&S;fqk#t17OX+I>2Ymis4W}qxty+>Wa^{eD5KoI*n`@QD"
        "tFLTLa$-f%LYt1t%_cQYL}->Xrn4lp#5d>Ru!G#x0llE%W2g&6F^)=BTDqOF`23=EvC^fO?SPOt67p6{3mSpZNU"
    ),
    "extract_all_hscan_regional_attenuation": (
        "c-qxkZFAeW5&o`Uf#Y>1WR9YpTrYFgl<DYVJNESAWE|x>b2^m+OQ6M?A{mgfBdh9f?_vRv07=Pidi~H$W-"
        "J1`i^bx>zLDR4bF?n`(RGp?vFwf%n`)V7-wzze8N6i~<FsPrd|1wDM&6AlBu-LRk~rst-7&r)$%@`Ea-UR7LJ4PdL9bIruFjB<&}"
        ">0&lSRfh<m&xJFc_S$n6A^RBr96YmufDeY#zNle);P7`=Hp2280aB<n0>}x61EmO0F5K2g0g~Woud`c}C_r=WH&1BDGq>reJZL%o"
        "CPX<T&_&EPgt2-"
        "<>~@`N~6<EYF6Ivi;a$drQ}4nSe1&rivDflQOwkQDFL>y&~riXu(0u^A(ArcwMbIiy|;rk@E`dlI0chmV<#F<2ME6CDYHGN$!GH%"
        "QU&xtRLX9=E>Hp0`xCQR_LJuYz_nPzgS2C<vb}iL6U0;X-"
        "*elexH=|j=>M5<vh>g<VF*{!6!Wf)>q><tcrLX!6MO$RdC!Dl|)N~NKGY?6XB>N=$y4CcEPZO(4t5;5ob4G6q-"
        "hcd86jAeI!MpQMq2NDBt)4ZvePDJsW@ec@kZX{}!EGT#*p%<bve^hEe1NoRxWc$6PM}Bg5gUKjerwVh@h^aYgSV@es{ZMl-"
        "f>YO(J>O+I{@D9KXar^!cBy@q}`|0wjQX&M2#lSOoO7D?m49f<(=h<mn(#=-"
        "__DmpHHhyb(A<BtyPTq`)9OinL9jVBlHFQdu$!|7G@@#5{fadiGGD2$LC3}Bhd@w?OL>J-"
        "K}3+DN%0PW>=xbyc9?q5FqxkwhTe_MR#MdINuDU&LyHU)b<9S?t{!-wPHpJvu1An=-"
        "T$>JzU7s4kX4SENQhR44JoZkPzkw*z?T)e$}e|7q1{P9#C<$}e4%$P+fy=JKkCR(!*i2=10`6nqW?vdeZ@pq&oz}3pcvm_NqQKWR"
        "vTqkgR;uMZC@*R$(>(U%N5`&SW*TzdIFs||K!f$a2M|kGA`XZTEQ&`d``Sm{lnlq`3bZQ}KEFIvGmJnog#XKZUV(>tf<dDq9X14W"
        "!_~*c~=#MT=oDtk~kfJU^P4SC33o?O&0t(0JoxQ|U*e3zcY(LFz(5s%UT9Q>RZYa3n`MPi&0F&c+jTBo?i?o3n9L>4_f9iG_YxO!"
        "SnP9Y+Y;CN=X+pH-+$k)?3HUiUJdAoS7ycpP1`g7stft~=W`hyj8e3@kfZwF~wd*|Rv4B1-"
        "qNLB4c@9uhH=5@OuH#CL?1jV+K>uZ)F|`s<0||hwCPx?c{-ShDfGN5lDa9Yx>tIcY9BA~)je=0U1I4X2AkZw-"
        "7t~5?%`mOy3ZT~yB2f?nVE0CTJqtQIQjQ}C@QRcbTLo!;F93k~gHx~+j@QBM^iiR>nY0R6Tm^Md2^TgsuyhuSH}w4q%(N*R!j1u("
        "sp)Y|Z?4qgFef#SbsC5c5vDaZR<yju!)W}F4_OBn(YRbf_Aw&Wx=7g+;Oz^H&k!}@7B490bYp-"
        "9JQ^5BB>}ZiX+;^260iWCse8*dK9Oz}4cR;kb7^q_GQe>u%mEcxk@Bh(d7<k>;9O@CV+Ch1gJ{06AY-"
        "m$LC0Rt&~6`u#L^acgyq^Z24KW=1krerz(x?rjfSv(zlw)dLzNy`w)Xzk>J3o3YE3P7-"
        "D|L?8%?WF1hhnW=+FrZuE^_<yhP)5=gZ`Qk=Nw7mh`mkaGI2i{32wWa-MTHb}n-"
        "yAEwDI108vSbwix5VO`LAt8uv#aao%L>|p5BS`Kc~^u-Kx{JduDfw^^@6?AfE-"
        "@wieYW3H;R*~x^C2s8HiLKL!;H4v56pg>zdg32#DXoO$0xis|i>Aiui6RJX0vi@eb})Hd@`0_|py}+pyFIeVS=pp4SmANVX2{lDa"
        "p%;6tu#lRS|(4-`%e#mt%_>X!Gz&61)<IxMRfvH*Vp+MR~9rQRm%8CH!^%8x=1;x-"
        "Q#`fL$#BZ5@`2J#+g(a=cT~l3X1@dNC<{tm1HjD4Wt>MUtoZIM_$N5A<&|d12vVk8M%hqPQI>+HKa9Q=Ll{#WFjhAHxp9d497Grm"
        "Q-cSa>A$?x|3^2Kd^Wmk?TB9<sjCv(ek0Dv=i=f&N8^KIT#^!B<<G7*EKg}BQMx1Z^vFa0tNL-gr#ukPL|6@$J~4a%h=iKoLR)<o"
        "h&4y>SnHsBPxp{b7^y%o&m+!i870@Z4++?TTt8DVUjbD5LT$dN)udP#J)hnmt;4Q2%RFId2(G~gj~(%HN`Am>MDSH<T^OuSV>lW?"
        "-hwL?9}&672u_jB<{Inr)n1pUxH_ql%*I!AtQ0~E9<H?DUH6C_A0ibaCedta%YG-E=sf@Zk6wH-Zur*g&5E0tSoi@y%+jcS}upw#"
        "c7Rx-5oJI5eV6+b3yb~>X4e(3Q)`MIhaHb9zsl$)Q5V@W=BDnrXpY2SFUQCr=ZUri%g02@1UNpa$T-"
        "q8rM)j;m}Kj*SIFD(@RXX!$(US6)dl*&^^Bgj;5%5{#<XxDC6D$B{@k%A1LP=QOMOAs)md3lAXHG<Pa@_t~tzF5rWzosHtg{9Ss5"
        "_Jeh<#hsH`l8$tRxr^7uK2m?+Na0U_FE1RwjPJ{MS$;J)jCg@K`{YW&p+AJRhoX;PCb3fD#_Q~5fSJDRq$Aqt0T{7!hNn4=ni&X#"
        "M`Ez~nzUqdF-m+|65C=+!YF&z_@-`_7@bwTZ8YBr4COjF^t7X-NSLkc5KSZH@LDcXje%p_sp~HTwx6dQoe?WE63m-"
        "WXZRYT?6R5qhvEcxtQmHB)k#04v=5S}pI0w5MWS)sYxUrhoBJG}*RHxm9YL3ir>98lf_OGeVTQ#dgmA|QEAAN>6)GApqP;6B|J_k"
        "92-_ZB;uN$wV>H322-g4Z>Wo4Q?9le}+wt_t=X*=3sjv|t206J(ZI>W$-{{2X(-"
        "zKPg&<MvOSTUNpEdf5EV9fBPXI9E8ft92bKtgZ8+BF4wtd^YV$R{J;;+z>sl+>z0k5&DX`uIko*DUAvWfRSyvfcmas_wsd;ZwX~@"
        "BKb>>(#N#Tqs`SiQ8{Y+?H)^a7spDA+&lU2g+m7MRhSX>C#UBWf*rRbG3|SGP9Q{SuM^=(G-"
        "6lc1>TTg7NBBbQ|SE8Df+mc)O$g79XW%>J3x_2SZe2+*UhWCdeaeTZV(M?AvjBeJ<~b?#^&^;?+sY;84u)-"
        "sj=8ciY%~Wxmd$YkV{C|CGU|$ttO^g|1rBhZ5Uj5uAVrI};!tnY$Djt!|9#mu9UkTOoWe)l)#GYX;q>?rrGTK@~Tsn?9}7rIDe3P"
        "UYJKWYYc370r?uLQiidB9Uvn_yb8shOBvmd3h%b(3mFa8n*v6(loD>rkbxp-"
        "`F|BF55O`nXC8W?eshQlE5UC>04zST>)z5_p)YO=)@X|#YU|*@worZda1MDaBG5g=NoV3!)+VWw(#D3mMF}cwg2{6o;}kohz`DG("
        "0i80@w1)9CL0;y_Q2)2*^%=H&$0YXevjQ+86!nE768Y-_&U=ReJl{SLO)+y-"
        "zV5E4__F~F>S)>7A1_1|G3y4*=L4hGom{s+CJi6#f=vgN-=D1_mtEYnG_Xv662<wfT#%7-"
        "!4G|hOqhBC3`~j{)#C6J4N%aj_Q^nv|9}8Up7&5TI|c56%{$fKLW<I6k{VozJ=X?25<h;&!@-"
        "7FT|HFj1{C>GJ@fQKph4dTwuw4dUL=A$F{SDzDj(5MQ-"
        ";g`pfdi%5zJE*TeT%iGH9|Z&w~quF1cqHvjKx^xvsgf2?NTTl4jrs!4xDY_W2mP7Wjv2N2vknBBlYTp2*Az%ucY`M!apnKocM-"
        "3zKbf*9snq1BguH9Vu>zCAnk3fy<3Sn(8lsgihF|IquJ1?;n{x9H7+d0uQ>>up9ITJhCH<7ew}rLh&h`!tq5-"
        "uW&22`{tGDg!If^(csM<DfEFLnw6@JTEe;<cg)i_TUp-&kn`!lo)L88kRIVyodHy=MHKf#8}(wypxejSl`r8ej~aRkR0l7-"
        "GWOvpt7KIm1*IiNdIbrF~5OmnN=TfrU_UJR?c~X-ZHeqm~;5EMNwR1!C!(Kr4?C12|}oRO<I4q=-"
        "3ilpE0&ya(RaC0bMM_Mp0?k84e*(4j}?MKB$b}K8E6Mw7ysK4vLU*^njG+^R&Z)<Q|6Tsz8*e6GK@5)MAf=9}bFW60}Fe{J39Cle"
        "!Iuh8r8A8ylh<8!jIV_hPJI06?zv_!n^gBZp+FljU8F3`Hn*FG9G9IcHG?j?|c3va|@DlcdB~*&;^QAB)DEeO|L{zEPiwg8N^P#a"
        "i6cJ)gA(j{i8eu=LO;fS-"
        "P7tBB&|Fs5`)0_;7zvM#d{VP0am*WfTB*;S3i8*3+m0w&3#$dOlB(A61;VnUgA>1*bc_(`Q84Rm4V4|>wdhO8O%_Vk62?1{}Gr2V"
        "=Mk>LFLHj(hVYVku!0iL#83x68MsZ&j97CYbc2LA@G#%Y-"
    ),
}


def _decode(payload: str) -> str:
    return zlib.decompress(base64.b85decode(payload.encode("ascii"))).decode("utf-8")


def _load_bundled_module(name: str, payload: str) -> None:
    module = types.ModuleType(name)
    module.__file__ = str(Path(__file__).resolve())
    module.__package__ = ""
    sys.modules[name] = module
    exec(compile(_decode(payload), module.__file__, "exec"), module.__dict__)


def _auto_sigma_hz(config) -> float:
    f_lo = (1.0 - config.bw_fc) * config.fc
    f_hi = (1.0 + config.bw_fc) * config.fc
    return (f_hi - f_lo) / 3.0


def _parser(default_raw: Path, default_output: Path, default_sigma: float) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate kidney-masked H-scan MAT files from RF data and kidney ROIs."
    )
    parser.add_argument("--raw-dir", type=Path, default=default_raw)
    parser.add_argument("--output-dir", type=Path, default=default_output)
    parser.add_argument("--kidney-alpha", type=float, default=0.5)
    parser.add_argument("--non-kidney-alpha", type=float, default=0.8)
    parser.add_argument(
        "--attenuation-taper-sigma-hz",
        type=float,
        default=default_sigma,
        help=f"Gaussian taper sigma in Hz (default: auto-derived {default_sigma:g}).",
    )
    parser.add_argument("--no-attenuation-taper", action="store_true")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    os.environ.setdefault("ANALYSIS_DIR", str(script_dir))
    for module_name, payload in _BUNDLED_MODULES.items():
        _load_bundled_module(module_name, payload)

    from config import Config
    import extract_all_hscan_regional_attenuation as extraction

    sigma_hz = _auto_sigma_hz(Config())
    project_dir = script_dir.parent
    default_raw = project_dir / "Data" / "Raw_data_cleaned"
    default_output = script_dir / (
        "HScan_all_valid_RF_kidneyROI_TGCremoved_Atten_"
        "kidney_0p5_non_kidney_0p8_sigma_auto1e6"
    )
    args = _parser(default_raw, default_output, sigma_hz).parse_args()

    print("H-scan settings:")
    print(f"  raw directory:       {args.raw_dir.resolve()}")
    print(f"  output directory:    {args.output_dir.resolve()}")
    print(f"  kidney attenuation:  {args.kidney_alpha} dB/(MHz cm)")
    print(f"  other attenuation:   {args.non_kidney_alpha} dB/(MHz cm)")
    print(f"  attenuation taper:   {not args.no_attenuation_taper}")
    print(f"  taper sigma:          {args.attenuation_taper_sigma_hz:g} Hz")
    print("  TGC removal:          True")

    extraction.run_batch(
        raw_dir=args.raw_dir,
        output_dir=args.output_dir,
        kidney_alpha=args.kidney_alpha,
        non_kidney_alpha=args.non_kidney_alpha,
        attenuation_taper_sigma_hz=args.attenuation_taper_sigma_hz,
        use_attenuation_taper=not args.no_attenuation_taper,
        limit=args.limit,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()

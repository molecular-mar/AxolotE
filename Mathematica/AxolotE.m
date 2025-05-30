(* ::Package:: *)

(* ::Section:: *)
(*Recursion relations for the expansion of RSHGTF product in HGTFs*)


(* ::Text:: *)
(*For the generation of the E coefficients, recurrence relations are used: one for increasing \[ScriptL] and another for incresing both \[ScriptL] and m of one of the RSHGTF. *)
(*All this relations are documented in: *)
(*(1) Pisani, C.; Dovesi, R.; Roetti, C. Hartree-Fock Ab Initio Treatment of Crystalline Systems, Lecture Notes in Chemistry, Vol. 48,; Lecture Notes in Chemistry; Springer Berlin Heidelberg, 1988; Vol. 48. *)
(*Note: X_0^-0 should be taken as 0; under certain circumstances, not considering this introduces problems.*)


(* ::Text:: *)
(*IMPORTANT NOTE: The parameters Dr, Dx, Dy, Dz, their "p"rime equivalents, and p, should be changed accordingly with the variable names used in the code generation.*)


(* ::Subsubsection:: *)
(*Recurrence  relationships  condensed  function*)


(* ::Text:: *)
(*Condensed function with the three recurrence relationships. It uses memoization, in order to used already produced and simplified coefficients as it runs.*)


(* ::Input:: *)
(*ExpanCll[l_,m_,lp_,mp_,t_,u_,v_]:=ExpanCll[l,m,lp,mp,t,u,v]=FullSimplify[*)
(*If[t<0 || u<0 || v<0 || t+u+v>l+lp||l<Abs[m]||lp<Abs[mp],0,*)
(*							If[t==0&& u==0 &&v==0 && l==0 &&lp==0,1,*)
(*							If[l!= 0 && lp ==0 ,ExpanC[l,m,t,u,v],*)
(*							If[lp != 0 && l ==0, ExpanCp[lp,mp,t,u,v],*)
(*							If[l>= lp,If[l==Abs[m] ,If[ l>1,If[m>0,*)
(*	(2*(l-1)+1)*(ExpanCll[l-1,l-1,lp,mp,t-1,u,v]/(2*p)-ExpanCll[l-1,-(l-1),lp,mp,t,u-1,v]/(2*p)+Dx*ExpanCll[l-1,l-1,lp,mp,t,u,v]-Dy*ExpanCll[l-1,-(l-1),lp,mp,t,u,v]+(t+1)*ExpanCll[l-1,l-1,lp,mp,t+1,u,v]-(u+1)*ExpanCll[l-1,-(l-1),lp,mp,t,u+1,v])*)
(*							,*)
(*	(2*(l-1)+1)*(ExpanCll[l-1,-(l-1),lp,mp,t-1,u,v]/(2*p)+ExpanCll[l-1,(l-1),lp,mp,t,u-1,v]/(2*p)+Dx*ExpanCll[l-1,-(l-1),lp,mp,t,u,v]+Dy*ExpanCll[l-1,(l-1),lp,mp,t,u,v]+(t+1)*ExpanCll[l-1,-(l-1),lp,mp,t+1,u,v]+(u+1)*ExpanCll[l-1,(l-1),lp,mp,t,u+1,v])*)
(*							],If[m>0,*)
(*	(2*(l-1)+1)*(ExpanCll[l-1,l-1,lp,mp,t-1,u,v]/(2*p)+Dx*ExpanCll[l-1,l-1,lp,mp,t,u,v]+(t+1)*ExpanCll[l-1,l-1,lp,mp,t+1,u,v])*)
(*							,*)
(*	(2*(l-1)+1)*(ExpanCll[l-1,(l-1),lp,mp,t,u-1,v]/(2*p)+Dy*ExpanCll[l-1,(l-1),lp,mp,t,u,v]+(u+1)*ExpanCll[l-1,(l-1),lp,mp,t,u+1,v])]],*)
(*	(1/(l-1-Abs[m]+1))*((2*(l-1)+1)*((ExpanCll[l-1,m,lp,mp,t,u,v-1]/(2*p))+(Dz*ExpanCll[l-1,m,lp,mp,t,u,v])+((v+1)*ExpanCll[l-1,m,lp,mp,t,u,v+1]))-((l-1+Abs[m])*((ExpanCll[l-2,m,lp,mp,t-2,u,v]+ExpanCll[l-2,m,lp,mp,t,u-2,v]+ExpanCll[l-2,m,lp,mp,t,u,v-2])/((2*p)^2)+((Dx*ExpanCll[l-2,m,lp,mp,t-1,u,v]+Dy*ExpanCll[l-2,m,lp,mp,t,u-1,v]+Dz*ExpanCll[l-2,m,lp,mp,t,u,v-1])/p)+((Dr2+((2*t+2*u+2*v+3)/(2*p)))*ExpanCll[l-2,m,lp,mp,t,u,v])+(2*Dx*(t+1)*ExpanCll[l-2,m,lp,mp,t+1,u,v])+(2*Dy*(u+1)*ExpanCll[l-2,m,lp,mp,t,u+1,v])+(2*Dz*(v+1)*ExpanCll[l-2,m,lp,mp,t,u,v+1])+*)
(*((t+2)*(t+1)*ExpanCll[l-2,m,lp,mp,t+2,u,v])+((u+2)*(u+1)*ExpanCll[l-2,m,lp,mp,t,u+2,v])+((v+2)*(v+1)*ExpanCll[l-2,m,lp,mp,t,u,v+2]))))*)
(*							]*)
(*							,*)
(*							If[lp==Abs[mp] ,If[ lp>1,If[mp>0,*)
(*	(2*(lp-1)+1)*(ExpanCll[l,m,lp-1,lp-1,t-1,u,v]/(2*p)-ExpanCll[l,m,lp-1,-(lp-1),t,u-1,v]/(2*p)+Dxp*ExpanCll[l,m,lp-1,lp-1,t,u,v]-Dyp*ExpanCll[l,m,lp-1,-(lp-1),t,u,v]+(t+1)*ExpanCll[l,m,lp-1,lp-1,t+1,u,v]-(u+1)*ExpanCll[l,m,lp-1,-(lp-1),t,u+1,v])*)
(*							,*)
(*	(2*(lp-1)+1)*(ExpanCll[l,m,lp-1,-(lp-1),t-1,u,v]/(2*p)+ExpanCll[l,m,lp-1,(lp-1),t,u-1,v]/(2*p)+Dxp*ExpanCll[l,m,lp-1,-(lp-1),t,u,v]+Dyp*ExpanCll[l,m,lp-1,(lp-1),t,u,v]+(t+1)*ExpanCll[l,m,lp-1,-(lp-1),t+1,u,v]+(u+1)*ExpanCll[l,m,lp-1,(lp-1),t,u+1,v])*)
(*							],If[mp>0,*)
(*	(2*(lp-1)+1)*(ExpanCll[l,m,lp-1,lp-1,t-1,u,v]/(2*p)+Dxp*ExpanCll[l,m,lp-1,lp-1,t,u,v]+(t+1)*ExpanCll[l,m,lp-1,lp-1,t+1,u,v])*)
(*							,*)
(*	(2*(lp-1)+1)*(ExpanCll[l,m,lp-1,(lp-1),t,u-1,v]/(2*p)+Dyp*ExpanCll[l,m,lp-1,(lp-1),t,u,v]+(u+1)*ExpanCll[l,m,lp-1,(lp-1),t,u+1,v])]],*)
(*	(1/(lp-1-Abs[mp]+1))*((2*(lp-1)+1)*((ExpanCll[l,m,lp-1,mp,t,u,v-1]/(2*p))+(Dzp*ExpanCll[l,m,lp-1,mp,t,u,v])+((v+1)*ExpanCll[l,m,lp-1,mp,t,u,v+1]))-((lp-1+Abs[mp])*((ExpanCll[l,m,lp-2,mp,t-2,u,v]+ExpanCll[l,m,lp-2,mp,t,u-2,v]+ExpanCll[l,m,lp-2,mp,t,u,v-2])/((2*p)^2)+((Dxp*ExpanCll[l,m,lp-2,mp,t-1,u,v]+Dyp*ExpanCll[l,m,lp-2,mp,t,u-1,v]+Dzp*ExpanCll[l,m,lp-2,mp,t,u,v-1])/p)+((Drp2+((2*t+2*u+2*v+3)/(2*p)))*ExpanCll[l,m,lp-2,mp,t,u,v])+(2*Dxp*(t+1)*ExpanCll[l,m,lp-2,mp,t+1,u,v])+(2*Dyp*(u+1)*ExpanCll[l,m,lp-2,mp,t,u+1,v])+(2*Dzp*(v+1)*ExpanCll[l,m,lp-2,mp,t,u,v+1])+*)
(*((t+2)*(t+1)*ExpanCll[l,m,lp-2,mp,t+2,u,v])+((u+2)*(u+1)*ExpanCll[l,m,lp-2,mp,t,u+2,v])+((v+2)*(v+1)*ExpanCll[l,m,lp-2,mp,t,u,v+2]))))*)
(*							]*)
(*							]*)
(*							]*)
(*							]*)
(*							]*)
(*]*)
(*]*)


(* ::Text:: *)
(*Generation of a Table containing the coefficients from l=l'=0 to l,l'=lmax,lpmax. It includes the cases where t+u+v  > l+l'. Used later for the actual code generation.*)


(* ::Input:: *)
(*SetOptions[$FrontEndSession,EvaluationCompletionAction->"ShowTiming"]*)


(* ::Input:: *)
(*lmax=3;*)
(*lpmax=3;*)
(*GenCoefs =Table[ExpanCll[li,mi,lpi,mpi,t,u,v],{li,0,lmax},{mi,-li,li},{lpi,0,lpmax},{mpi,-lpi,lpi},{t,0,li+lpi},{u,0,li+lpi-t},{v,0,li+lpi-t-u}];*)


(* ::Subsubsection:: *)
(*Following Section has different tests which are not necessary to run for actual code generation.*)


(* ::Text:: *)
(*Note on GenCoefs: When accessing a given table element, you need to adapt the desired indices (l,m,l',p',t,u,v) to the actual indices of*)
(*the table. For that purpose you can use the following transformations (X id correspond to the table index):*)
(*l id-> l+1 ------------------------------- l' id-> l'+1*)
(*m id ->m+(l+1)------------------------m' id ->m'+(l'+1)*)
(*t|u|v id -> t+1|u+1|v+1*)
(**)
(*For example take the following commands and the relations used for accessing GenCoefs elements.*)


(* ::Text:: *)
(*To print a section of the table for a given set of lm,l'm':*)


(* ::Input:: *)
(*(*Adjust Print variables to print the desired combination. Check that the value was previously calculated.*)*)
(*lPrint=1;*)
(*mPrint=0;*)
(*lpPrint=0;*)
(*mpPrint=0;*)
(*FlattCoeff=Flatten[GenCoefs[[lPrint+1,mPrint+(lPrint+1),lpPrint+1,mpPrint+(lpPrint+1)]]];*)
(*FlattTUV=Flatten[Table[ToString[t]<>ToString[u]<>ToString[v],*)
(*{t,0,lPrint+lpPrint},{u,0,lPrint+lpPrint-t},{v,0,lPrint+lpPrint-t-u}]];*)
(*Print["Table for lm,l'm': ",lPrint,mPrint,",",lpPrint,mpPrint]*)
(*Grid[{FlattTUV,FlattCoeff},Frame->All]*)


(* ::Text:: *)
(*To print all tables up to a given set of lm,l'm' values:*)


(* ::Input:: *)
(*Nvals =0;(*To store the number of coefficients*)*)
(*maxL=1;*)
(*maxLp=1;*)
(*For[l=0,l<=maxL ,l+=1,*)
(*For[lp=0,lp<=maxLp,lp+=1,*)
(*For[m=-l,m<=  l,m+=1,	*)
(*For[mp=-lp,mp<=  lp,mp+=1,*)
(**)
(*FlattCoeff=Flatten[GenCoefs[[l+1,m+(l+1),lp+1,mp+(lp+1)]]];*)
(*FlattTUV=Flatten[Table[ToString[t]<>ToString[u]<>ToString[v],*)
(*{t,0,l+lp},{u,0,l+lp-t},{v,0,l+lp-t-u}]];*)
(*Print["Table for lm,l'm': ",l,m,",",lp,mp];*)
(*Print[Grid[{FlattTUV,FlattCoeff},Frame->All]];*)
(*Nvals=Nvals+1*)
(*	]*)
(*	]*)
(*	]*)
(*	]*)
(*	*)


(* ::Input:: *)
(*(*This code section might be used to recover the common factors for a given l,l' combination*)*)
(*(*MemberQ[GenCoefs,p,Infinity]*)
(*MemberQ[{pDx*pDy,1/p,pDz},pDz*_] || MemberQ[{pDx*pDy,p,pDz},pDz]*)
(*MemberQ[{pDx*pDy,pDz},p,Infinity]*)
(*ContainsAny[{pDx*pDy,1/p,pDz},{pDz*_,pDz}]*)
(*!FreeQ[{pDx*pDy,1/p,pDz},p]*)
(*MemberQ[{pDx*pDy,1/(2*p),pDz},p^_*_]*)
(*(*The result below could be used to recover all the unique values in the table.*)*)
(*DeleteDuplicates@*Flatten@GenCoefs*)*)
(**)


(* ::Text:: *)
(*The parameters should be changed accordingly to the variables used inside the generated code. Example:*)


(* ::Input:: *)
(*cVal=CForm[FullSimplify[GenCoefs[[2,1,2,1,1,3,1]]]];*)
(*GenCoefs[[2,1,2,1,1,3,1]]*)
(*GenCoefs[[2,1,2,1]]*)


(* ::Input:: *)
(*(*To change Power to std::pow*)*)
(*StringReplace[ToString[cVal],"Power"->"std::pow"]*)


(* ::Input:: *)
(*(*Replace expressions. NOTE: Use before changing the definition of Dx,Dy,Dz,Dr2. This happens below.*)*)


(* ::Input:: *)
(*GenCoefs[[2,1,2,1,1,2,1]]/.{Dyp->yPToB ,Dy->yPToA}*)


(* ::Subsubsection:: *)
(*Autogeneration  of  C++  code*)


(* ::Text:: *)
(*Given the table generated, print C++ functions for the necessary indexes. Works for a given pair of l and l' values each time. *)


(* ::Input:: *)
(*maxL=4;*)
(*maxLp=4;*)
(*Nvals =0;*)
(*(*Use OpenWrite to start a new file (previous content will be erased.*)*)
(*(*IMPORTANT: write here a valid path with the desired final filename.*)*)
(*ECoefFile = OpenWrite["/home/user/AxolotE_Mat.cpp"]*)
(*(*Use OpenAppend to write to an existing file. Remember to specify correct file name*)
(*ECoefFile = OpenAppend["/home/user/AxolotE_Mat.cpp"]*)*)
(*(*Function title, with generic name ECoeffLAndL'. New line for avoiding overlap of functions.*)*)
(*BaseName="ECoeff";*)
(*(*This line is added to stablish the usage of the routines inside OpenACC parallel regions*)*)
(*OpenACC="\n#pragma acc routine seq";*)
(**)
(*For[l=0,l<=maxL ,l+=1,*)
(*For[lp=0,lp<=maxLp,lp+=1,*)
(*Nvals =0;*)
(*WriteString[ECoefFile,OpenACC,"\nstd::vector<double> ",BaseName,l,"And",lp,"(double expAlpha, double expBeta,\[IndentingNewLine]\tstd::array<double,3> muCoords, std::array<double,3> nuCoords)\[IndentingNewLine]{\[IndentingNewLine]\tdouble expGamma{expAlpha + expBeta};*)
(*\tdouble betaOverGamma{expBeta / expGamma};*)
(*\tdouble xPtoA {betaOverGamma*(nuCoords[0] - muCoords[0])};*)
(*\tdouble yPtoA {betaOverGamma*(nuCoords[1] - muCoords[1])};*)
(*\tdouble zPtoA {betaOverGamma*(nuCoords[2] - muCoords[2])};*)
(*\tdouble d2PtoA {xPtoA*xPtoA + yPtoA*yPtoA + zPtoA*zPtoA};*)
(*\tdouble alphaOverGamma{-expAlpha / expGamma};*)
(*\tdouble xPtoB {alphaOverGamma*(nuCoords[0] - muCoords[0])};*)
(*\tdouble yPtoB {alphaOverGamma*(nuCoords[1] - muCoords[1])};*)
(*\tdouble zPtoB {alphaOverGamma*(nuCoords[2] - muCoords[2])};*)
(*\tdouble d2PtoB {xPtoB*xPtoB + yPtoB*yPtoB + zPtoB*zPtoB};*)
(*\tint nCoeffs {",(2*l+1)*(2*lp+1)*(l+lp+1)*(l+lp+2)*(l+lp+3)/6,"};*)
(*\tstd::vector<double> expCoeffs(nCoeffs,0);"]*)
(*For[m=-l,m<=  l,m+=1,*)
(*For[mp=-lp,mp<=  lp,mp+=1,*)
(*For[v=0,v<=maxL+maxLp,v++,*)
(*	For[u=0,u<=maxL+maxLp,u++,*)
(*	For[t=0,t<=maxL+maxLp,t++,*)
(*If[t+u+v<=l+lp,*)
(*(*Recover the coefficient from the table*)*)
(*currentCoef =CForm[FullSimplify[GenCoefs[[l+1,m+(l+1),lp+1,mp+(lp+1),v+1,u+1,t+1]]]/.\*)
(*{Dx->xPtoA, Dy-> yPtoA,Dz->zPtoA,\*)
(*Dxp->xPtoB, Dyp-> yPtoB,Dzp->zPtoB,\*)
(*Dr2-> d2PtoA, Drp2 -> d2PtoB, p->expGamma}];*)
(*WriteString[ECoefFile,"\n\t//"," m,mp=",m,mp," t,u,v=",t,u,v];	*)
(*WriteString[ECoefFile,"\n\texpCoeffs[",Nvals,"] = ",*)
(*StringReplace[ToString[currentCoef],"Power"->"std::pow"] ,";"];*)
(*Nvals=Nvals+1*)
(*	]]*)
(*	]*)
(*	]*)
(*]*)
(*]*)
(*(*Return and closing section*)*)
(*WriteString[ECoefFile,"\n\treturn expCoeffs;\[IndentingNewLine]}\n"];*)
(*]*)
(*]*)
(**)
(*Close[ECoefFile]*)


(* ::Section:: *)
(*Testing the generated coefficients*)


(* ::Text:: *)
(*We need to test if the coefficients are correct. For this, the explicit relation of the linear combination is evaluated using the generated expressions.*)


(* ::Input:: *)
(*Function to transform expresion in two spherical coordinates sets to cartesian . Note the variables used for the coordinates: rA and rB .*)


(* ::Input:: *)
(*TtoCart[SpherFun_]:=TransformedField["Spherical"->"Cartesian",*)
(*TransformedField["Spherical"->"Cartesian",*)
(*FullSimplify[SpherFun,*)
(*{rA>= 0 &&rB>=  0 *)
(*&&0<=\[Theta]A<=\[Pi] &&0<=\[Theta]B<=\[Pi] *)
(*&& 0<=\[Phi]A<2\[Pi] && 0<= \[Phi]B<2\[Pi]}],*)
(*{rA,\[Theta]A,\[Phi]A}->{xA,yA,zA}],*)
(*{rB,\[Theta]B,\[Phi]B}->{xB,yB,zB}*)
(*]*)


(* ::Text:: *)
(*Real Solid Harmonics *)


(* ::Input:: *)
(*rDist[x_,y_,z_]:=Sqrt[x^2+y^2+z^2];*)
(*rSphericalH[r_,l_,m_,\[Theta]_,\[Phi]_]:=r^l*(1/((-1)^m))*LegendreP[l,Abs[m],Cos[\[Theta]]]*Exp[I*m*\[Phi]];*)
(*rRealSphericalH[r_,l_,m_,\[Theta]_,\[Phi]_]:=Piecewise[{{(rSphericalH[r,l,m,\[Theta],\[Phi]]+rSphericalH[r,l,-m,\[Theta],\[Phi]])/(2),m>= 0},{(rSphericalH[r,l,-m,\[Theta],\[Phi]]-rSphericalH[r,l,m,\[Theta],\[Phi]])/(2*I),m< 0}}]*)
(*(*Analog definition*)*)
(*RSolidH2[L_,m_,\[Theta]_,\[Phi]_,r_]:=r^L*LegendreP[L,Abs[m],Cos[\[Theta]]]*Piecewise[{{Sin[Abs[m]*\[Phi]],m<0},{Sqrt[1/2],m==0},{Cos[Abs[m]*\[Phi]],m>0}}];*)


(* ::Text:: *)
(*Table of Real Solid Harmonics in spherical and cartesian coordinates. Different approaches for simplifying the expressions are shown.*)


(* ::Input:: *)
(*(*Checked for lmaxP=4 (associated with g orbitals)*)*)
(*lmaxP=4;*)
(*For[i=0,i<=lmaxP,i++,*)
(*Print[Style["-------------------------" ,18,Bold]];*)
(*Print[Style["L value: " ,Bold],i ];*)
(*Print[Style["-------------------------",18,Blue] ];*)
(*For[j=-i,j<=i,j++,*)
(*Print[Style["M value: ",Italic] ,j ];*)
(*(*Spherical*)*)
(*currentSH = rRealSphericalH[rA,i,j,\[Theta]A,\[Phi]A];*)
(*fcurrentSH = FullSimplify[currentSH];*)
(*(*Cartesian*)*)
(*tCurrentSH = TtoCart[rRealSphericalH[rA,i,j,\[Theta]A,\[Phi]A]];*)
(*tfcurrentSH = TtoCart[fcurrentSH];*)
(*fCurrentSH = FullSimplify[tCurrentSH];*)
(*FfCurrentSH = FullSimplify[tfcurrentSH];*)
(*(*Note that for l=4 m=4 only with TrigExpand the simplification gives the expected result*)*)
(*trigCurrentSH = FullSimplify[TrigExpand[tCurrentSH]];*)
(*Print[Style["Spherical",Italic]  ];*)
(*Print[currentSH];*)
(*Print[fcurrentSH];*)
(*Print[Style["Cartesian",Italic]  ];*)
(*Print[tCurrentSH];*)
(*Print[tfcurrentSH];*)
(*Print[fCurrentSH];*)
(*Print[FfCurrentSH];*)
(*Print[trigCurrentSH];*)
(*Print["***********"]*)
(*];*)
(*Print[Style["-------------------------" ,18,Blue]];*)
(*]*)


(* ::Text:: *)
(*Product of Real Solid Harmonics in spherical coordinates (not used later):*)


(* ::Input:: *)
(*productRSH [r1_,l1_,m1_,\[Theta]1_,\[Phi]1_,r2_,l2_,m2_,\[Theta]2_,\[Phi]2_]:=Assuming[r1>=  0 &&r2>=  0 &&0<=\[Theta]1<=\[Pi]/2 &&0<=\[Theta]2<=\[Pi]/2 && 0<=\[Phi]1<=\[Pi] && 0<=\[Phi]2<=\[Pi] ,rRealSphericalH[r1,l1,m1,\[Theta]1,\[Phi]1]*rRealSphericalH[r2,l2,m2,\[Theta]2,\[Phi]2]]*)


(* ::Text:: *)
(*Product of Real Solid Harmonics in cartesian coordinates:*)


(* ::Input:: *)
(*cartProdRSH[rA_,lA_,mA_,\[Theta]A_,\[Phi]A_,rB_,lB_,mB_,\[Theta]B_,\[Phi]B_]:= FullSimplify[TrigExpand[TtoCart[rRealSphericalH[rA,lA,mA,\[Theta]A,\[Phi]A]]]]*FullSimplify[TrigExpand[TtoCart[rRealSphericalH[rB,lB,mB,\[Theta]B,\[Phi]B]]]]*)


(* ::Text:: *)
(*Coordinates for overlap distribution:*)


(* ::Input:: *)
(*pCoor[\[Alpha]_,rA_,\[Beta]_,rB_]:=(\[Alpha]*rA+\[Beta]*rB)/(\[Alpha]+\[Beta]);*)
(*gamma[\[Alpha]_,\[Beta]_]:=\[Alpha]+\[Beta];*)
(*pComp[\[Alpha]_,\[Beta]_,compA_,compB_]:=(\[Alpha]*compA+\[Beta]*compB)/(\[Alpha]+\[Beta]);*)


(* ::Text:: *)
(*Distance from P to A or B (Nuclear positions, not to be confuse with xA coordinates).*)


(* ::Input:: *)
(*DToAorB[AorB_]:=pCoor[\[Alpha],A,\[Beta],B]-AorB/.{xA->x-Ax,xB-> x-Bx,*)
(*yA->y-Ay,yB-> y-By,*)
(*zA->z-Az,zB-> z-Bz};*)
(*DComp[compA_,compB_,compD_]:=pComp[\[Alpha],\[Beta],compA,compB]-compD/.{xA->x-Ax,xB-> x-Bx};*)


(* ::Text:: *)
(*Hermite polynomials, adjusted as they appear on Hermite Gaussians.*)


(* ::Input:: *)
(*singleHermite [\[Alpha]_,x_,id_]:=\[Alpha]^(id/2)*HermiteH[id,x*Sqrt[\[Alpha]]];*)
(*fullHermite[\[Alpha]_,x_,y_,z_,t_,u_,v_]:=singleHermite[\[Alpha],x,t]*singleHermite[\[Alpha],y,u]*singleHermite[\[Alpha],z,v];*)


(* ::Text:: *)
(*Using the coefficients relations obtained, sum over t+u+v <= l+l'.*)


(* ::Input:: *)
(*sumHermite[\[Gamma]_,x_,y_,z_,l1_,m1_,l2_,m2_]:=Sum[ExpanCll[l1,m1,l2,m2,t,u,v]*fullHermite[\[Gamma],x,y,z,t,u,v]Boole[t+u+v<=l1+l2],{t,0,l1+l2},{v,0,l1+l2},{u,0,l1+l2}];*)


(* ::Input:: *)
(*pExp=gamma[\[Alpha],\[Beta]];*)
(*pX=pComp[\[Alpha],\[Beta],x-Ax,x-Bx];*)
(*pY=pComp[\[Alpha],\[Beta],y-Ay,y-By];*)
(*pZ=pComp[\[Alpha],\[Beta],z-Az,z-Bz];*)
(*Dx=DComp[Ax,Bx,Ax];*)
(*Dy=DComp[Ay,By,Ay];*)
(*Dz=DComp[Az,Bz,Az];*)
(*Dr2 = (pComp[\[Alpha],\[Beta],Ax,Bx]-Ax)^2+(pComp[\[Alpha],\[Beta],Ay,By]-Ay)^2+(pComp[\[Alpha],\[Beta],Az,Bz]-Az)^2;*)
(*Dxp=DComp[Ax,Bx,Bx];*)
(*Dyp=DComp[Ay,By,By];*)
(*Dzp=DComp[Az,Bz,Bz];*)
(*Drp2 = (pComp[\[Alpha],\[Beta],Ax,Bx]-Bx)^2+(pComp[\[Alpha],\[Beta],Ay,By]-By)^2+(pComp[\[Alpha],\[Beta],Az,Bz]-Bz)^2; *)
(*p=pExp;*)


(* ::Input:: *)
(*Dr2*)
(*Drp2*)


(* ::Input:: *)
(*(*To check the expression for the lineal combination of Hermite polynomials with the calculated coefficients*)*)
(*Sum[fullHermite[pExp,xP,yP,zP,t,u,v]Boole[t+u+v<=4+4],{t,0,4+4},{v,0,4+4},{u,0,4+4}]*)
(*(*Number of terms, including those where the coefficient is 0*)*)
(*Sum[1Boole[t+u+v<=4+4],{t,0,4+4},{v,0,4+4},{u,0,4+4}]*)


(* ::Input:: *)
(*(*Checked for up to lMax=4*)*)


(* ::Text:: *)
(*Testing the generated expressions for the coefficients, multiplying each with its corresponding Hermite polynomial and comparing the result with the associated product of Real Solid Harmonics. In a single core (Core i7 2.2GHz) and testing each possible (l,m,l',m') combination this test takes around 47 hours for maximum l and l' value of 4 (g orbitals).    *)


(* ::Input:: *)
(*maxLVal=4;*)
(*minLVal = 0;*)
(*minLpVal=0;*)
(*For[l=minLVal,l<= maxLVal,l++,*)
(*For[lp=minLpVal,lp<= maxLVal,lp++,*)
(*For[m=-l,m<= l,m++,*)
(*For[mp=-lp,mp<= lp,mp++,*)
(*Print[l,m,lp,mp];*)
(*hermite1 = sumHermite[pExp,pX,pY,pZ,l,m,lp,mp];*)
(*(*FullSimplify for hermite1 can fail or be too slow. Expanding and Together works better.*)*)
(*hermite2 = Together[Expand[hermite1]];*)
(*HermLessHarm = hermite2 -FullSimplify[cartProdRSH[rA,l,m,\[Theta]A,\[Phi]A,rB,lp,mp,\[Theta]B,\[Phi]B]]/.{xA->(x-Ax),xB-> (x-Bx),yA->(y-Ay),yB-> (y-By),zA->(z-Az),zB-> (z-Bz)};*)
(*Print[Simplify[HermLessHarm]]*)
(*]*)
(*]*)
(*]*)
(*]*)


(* ::Input:: *)
(*(*Given the memoization technique, this needs to be done to clear the stored values if needed.*)*)
(*Clear[ExpanCll]*)

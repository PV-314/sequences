\\ \r sequences\primitive-divisors\Lehmer-t-5-prim-div-checks.gp

read("sequences\\primitive-divisors\\Lehmer-utils.gp");

\\ use \p 50 (or higher)
\\ this, and tol=10.0^(-30), will work with bnd up to 10^8

\\ the point of the code in this file is to do searches for "small" values of a
\\ to find all Lehmer sequences with u_5 not having a primitive divisor
\\ the means of doing so is Cam's cyclotomic criterion that we must have \Phi_5=\pm 1, \pm 5
\\ \Phi_5 is a quadratic binary form in a and b, so for each value of a, we can solve
\\ for b
\\ then we check that the pair (a,b) fits with BHV/my paper

\\ entry point. Call this function. From outside, the other functions can be ignored.
\\ 7 May 2024
t5_check(dbg=0)={
	my(bnd);
	
	bnd=10^7;
	for (a=-bnd,bnd,
		check_a_c(a,1,dbg);
		check_a_c(a,-1,dbg);
		check_a_c(a,5,dbg);
		check_a_c(a,-5,dbg);
	);
}

\\ the expressions for b come from f:=5*a*a+10*a*b+b*b:solve(f=16*c,b);
\\ f is the numerator here (and the 16 in the 16*c on the RHS is the denominator here):
\\ restart:with(numtheory):al:=(sqrt(a)+sqrt(b))/2:be:=(sqrt(a)-sqrt(b))/2:
\\ myPhi5:=expand(y^4*cyclotomic(5,x/y));simplify(subs(x=al,y=be,myPhi5));print(numer(%),denom(%)):
\\ 9 May 2024
check_a_c(a,c,dbg=0)={
	my(b);

	if(issquare(5*a*a+4*c),
		b=-5*a+2*sqrtint(5*a*a+4*c);
		if(dbg!=0 && abs(a)<4,
			print("in check_a_c(): a=",a,", b=",b);
		);
		check_alpha_beta(a,b,dbg);
		if(5*a*a+4*c!=0,
			b=-5*a-2*sqrtint(5*a*a+4*c);
			if(dbg!=0 && abs(a)<4,
				print("in check_a_c(): a=",a,", b=",b);
			);
			check_alpha_beta(a,b,dbg);
		);
	);
}

\\ make alpha and beta from a and b (as in BHV)
\\ 9 May 2024
check_alpha_beta(a,b,dbg=0)={
	my(al,be,phi5,tol);
	
	tol=10.0^(-30);
	if(dbg!=0,
		print("in check_alpha_beta(): checking a=",a,", b=",b);
	);
	al=(sqrt(a)+sqrt(b))/2;
	be=(sqrt(a)-sqrt(b))/2;
	if(dbg!=0,
		print("in check_alpha_beta(): al=",al);
		print("in check_alpha_beta(): be=",be);
	);
	if(abs(be)>tol && (abs(al/be-1)>tol && abs(al/be+1)>tol && abs(al/be-I)>tol && abs(al/be+I)>tol && abs((al/be)^12-1)>tol),
		phi5=(5*a*a+10*a*b+b*b)/16;
		if(dbg!=0,
			print("in check_alpha_beta(): checking Phi_5=",phi5," for a=",a,", b=",b);
		);
		if(abs(phi5)==1 || abs(phi5)==5,
			if(dbg!=0,
				print("in check_alpha_beta(): checking sequence for a=",a,", b=",b,", Phi_5=",phi5);
			);
			check_a_b(a,b,phi5,dbg);
		);
	);
}

\\ do a and b come from where I want/expect?
\\ 7 May 2024
check_a_b(aArg,bArg,phi5,dbg=0)={
	my(a,b,b1a,b1b,e,isFound,k,kMin,n,nArray,p,q,sgn);
	
	kMin=-2;

	a=aArg;
	b=bArg;
	sgn="+";
	if(aArg<0,
		a=-aArg;
		b=-bArg;
		sgn="-";
	);
	if(abs(phi5)==5 && aArg==-1 && bArg==-5,
		a=aArg;
		b=bArg;
		sgn="+";
	);
	
	p=aArg;
	if((p-bArg)%4!=0,
		print("BAD: a=",aArg,", b=",bArg,", p=",p,", q=",(p-bArg)/4);
		1/0;
	);
	if((p-bArg)%4==0,
		q=(p-bArg)/4;
		if(abs(phi5)==1,
			nArray=get_fibonacci_index(a);
			if(dbg!=0,
				print("in check_a_b(): a=",a,", aArg=",aArg,", b=",b,", bArg=",bArg,", Phi_5=",phi5,", nArray=",nArray);
			);
			isFound=0;
			for(i=1,length(nArray),
				n=nArray[i];
				e=1;
				k=n+2*e;
				if(k<kMin || k-2*e<kMin,
					printf("   WARNING: index<%2d Fibon not yet implemented: k=%2d, e=%2d, k-2e=%2d: a=%8d, b=%10d, p=%8d, q=%9d, a=F_n for n=%2d, Phi_5=%2d\n",kMin,k,e,k-2*e,a,b,p,q,n,phi5);
				);
				if(k>kMin-1 && k-2*e>kMin-1,
					b1a=fibonacci(k-2*e)-4*fibonacci(k);
					if(b1a==b,
						printf("FOUND Fibon: a=%8d, b=%10d, p=%8d, q=%9d, (a,b)=%s(phi_(k-2e), phi_(k-2e)-4phi_k) for k=%2d, e=%2d, Phi_5=%2d\n",aArg,bArg,p,q,sgn,k,e,phi5);
						isFound=1;
					);
				);
				e=-1;
				k=n+2*e;
				if(k<kMin || k-2*e<kMin,
					printf("   WARNING: index<%2d Fibon not yet implemented: k=%2d, e=%2d, k-2e=%2d: a=%8d, b=%10d, p=%8d, q=%9d, a=F_n for n=%2d, Phi_5=%2d\n",kMin,k,e,k-2*e,a,b,p,q,n,phi5);
				);
				if(k>kMin-1 && k-2*e>kMin-1,
					b1b=fibonacci(k-2*e)-4*fibonacci(k);
					if(b1b==b,
						printf("FOUND Fibon: a=%8d, b=%10d, p=%8d, q=%9d, (a,b)=%s(phi_(k-2e), phi_(k-2e)-4phi_k) for k=%2d, e=%2d, Phi_5=%2d\n",aArg,bArg,p,q,sgn,k,e,phi5);
						isFound=1;
					);
				);
			);
			if(isFound==0,
				printf("\n(*) NOT FOUND b in Fibon: a=%8d, b=%10d, p=%8d, q=%9d, a=F_n for n=%3d, b1a=%3d, b1b=%3d, Phi_5=%2d\n\n",a,b,p,q,n,b1a,b1b,phi5);
				1/0;
			);
			if(length(nArray)==0,
				printf("\n(*) NOT FOUND a in Fibon: a=%8d, b=%10d, p=%8d, q=%9d, Phi_5=%2d\n\n",a,b,p,q,phi5);
				1/0;
			);
		);
		if(abs(phi5)==5,
			nArray=get_lucas_index(a);
			isFound=0;
			for(i=1,length(nArray),
				n=nArray[i];
				e=1;
				k=n+2*e;
				if(k<kMin || k-2*e<kMin,
					printf("   WARNING: index<%2d Lucas not yet implemented: k=%2d, e=%2d, k-2e=%2d: a=%8d, b=%10d, p=%8d, q=%9d, a=L_n for n=%2d, Phi_5=%2d\n",kMin,k,e,k-2*e,a,b,p,q,n,phi5);
				);
				if(k>kMin-1 && k-2*e>kMin-1,
					b1a=get_lucas(k-2*e)-4*get_lucas(k);
					if(b1a==b,
						printf("FOUND Lucas: a=%8d, b=%10d, p=%8d, q=%9d, (a,b)=%s(psi_(k-2e), psi_(k-2e)-4psi_k) for k=%2d, e=%2d, Phi_5=%2d\n",aArg,bArg,p,q,sgn,k,e,phi5);
						isFound=1;
					);
				);
				e=-1;
				k=n+2*e;
				if(k<kMin || k-2*e<kMin,
					printf("   WARNING: index<%2d Lucas not yet implemented: k=%2d, e=%2d, k-2e=%2d: a=%8d, b=%10d, p=%8d, q=%9d, a=L_n for n=%2d, Phi_5=%2d\n",kMin,k,e,k-2*e,a,b,p,q,n,phi5);
				);
				if(k>kMin-1 && k-2*e>kMin-1,
					b1b=get_lucas(k-2*e)-4*get_lucas(k);
					if(b1b==b,
						printf("FOUND Lucas: a=%8d, b=%10d, p=%8d, q=%9d, (a,b)=%s(psi_(k-2e), psi_(k-2e)-4psi_k) for k=%2d, e=%2d, Phi_5=%2d\n",aArg,bArg,p,q,sgn,k,e,phi5);
						isFound=1;
					);
				);
			);
			if(isFound==0,
				printf("\n(*) NOT FOUND b in Lucas: a=%8d, b=%10d, p=%8d, q=%9d, a=L_n for n=%3d, b1a=%2d, b1b=%2d, Phi_5=%2d\n\n",a,b,p,q,n,b1a,b1b,phi5);
				1/0;
			);
			if(length(nArray)==0,
				printf("\n(*) NOT FOUND a in Lucas: a=%8d, b=%10d, p=%8d, q=%9d, Phi_5=%2d\n\n",a,b,p,q,phi5);
				1/0;
			);
		);
	);
}

\\ returns [] if not found
\\ only looks for non-negative indices
\\ 7 May 2024
get_fibonacci_index(f)={
	my(i);
	
	if(f==-1,return([-2]));
	if(f==1,return([-1,1,2]));
	i=get_index(f,0,1,1,1,"Fibon(phi)");
	if(i!=-1000,
		return([i]);
	);
	return([]);
}

\\ returns [] if not found
\\ finds negative indices too
\\ 7 May 2024
xxx_real_get_fibonacci_index(f)={
	for(k=0,100,
		if(f==fibonacci(k),
			return([k]);
		);
		if(f==(-1)^(k+1)*fibonacci(k),
			return([-k]);
		);
	);
	return([]);
}

\\ 9 May 2024
output_lucas_sequence()={
	output_sequence(2,1,1,1);
}

\\ 7 May 2024
get_lucas(k)={
	if(k==-2,return(3));
	if(k==-1,return(-1));
	return(get_kth_element(k,2,1,1,1,"Lucas(psi)"));
}

\\ returns [] if not found
\\ only looks for non-negative indices
\\ 7 May 2024
get_lucas_index(f)={
	my(i);
	
	if(f==-1,return([-1]));
	if(f==3,return([-2,2]));
	i=get_index(f,2,1,1,1,"Lucas");
	if(i!=-1000,
		return([i]);
	);
	return([]);
}
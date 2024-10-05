\\ \r sequences\primitive-divisors\Lehmer-t-8-prim-div-checks.gp

read("sequences\\primitive-divisors\\Lehmer-utils.gp");

\\ use \p 50 (or higher)
\\ this, and tol=10.0^(-30), will work with bnd up to 10^8

\\ the point of the code in this file is to do searches for "small" values of a
\\ to find all Lehmer sequences with u_8 not having a primitive divisor
\\ the means of doing so is Cam's cyclotomic criterion that we must have \Phi_8=\pm 1, \pm 2
\\ \Phi_8 is a quadratic binary form in a and b, so for each value of a, we can solve
\\ for b
\\ then we check that the pair (a,b) fits with BHV/my paper

\\ entry point. Call this function. From outside, the other functions can be ignored.
\\ 7 May 2024
t8_check(dbg=0)={
	my(bnd);
	
	bnd=10^7;
	for (a=-bnd,bnd,
		check_a_c(a,1,dbg);
		check_a_c(a,-1,dbg);
		check_a_c(a,2,dbg);
		check_a_c(a,-2,dbg);
	);
}

\\ the expressions for b come from phi8:=a*a+6*a*b+b*b:solve(phi8=8*c,b);
\\ f is the numerator here (and the 8 in the 8*c on the RHS is the denominator here):
\\ restart:with(numtheory):al:=(sqrt(a)+sqrt(b))/2:be:=(sqrt(a)-sqrt(b))/2:
\\ myPhi8:=expand(y^4*cyclotomic(8,x/y));simplify(subs(x=al,y=be,myPhi8));print(numer(%),denom(%)):
\\ 9 May 2024
check_a_c(a,c,dbg=0)={
	my(b);

	if(issquare(2*a*a+2*c),
		b=-3*a+2*sqrtint(2*a*a+2*c);
		check_alpha_beta(a,b,dbg);
		if(2*a*a+2*c!=0,
			b=-3*a-2*sqrtint(2*a*a+2*c);
			check_alpha_beta(a,b,dbg);
		);
	);
}

\\ make alpha and beta from a and b (as in BHV)
\\ 9 May 2024
check_alpha_beta(a,b,dbg=0)={
	my(al,be,phi8,tol);
	
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
		phi8=(a*a+6*a*b+b*b)/8;
		if(dbg!=0,
			print("in check_alpha_beta(): checking Phi_8=",phi8," for a=",a,", b=",b);
		);
		if(abs(phi8)==1 || abs(phi8)==2,
			if(dbg!=0,
				print("in check_alpha_beta(): checking sequence for a=",a,", b=",b,", Phi_8=",phi8);
			);
			check_a_b(a,b,phi8,dbg);
		);
	);
}

\\ determine p and q. Then find them in the sequences
\\ do they come from where I want/expect?
\\ 7 May 2024
check_a_b(aArg,bArg,phi8,dbg=0)={
	my(a,b,b1a,b1b,e,isFound,k,n,nArray,p,q,sgn);
	a=aArg;
	b=bArg;
	sgn="+";
	if(aArg<0,
		a=-aArg;
		b=-bArg;
		sgn="-";
	);
	
	p=aArg;
	if((p-bArg)%4!=0,
		print("BAD: a=",aArg,", b=",bArg,", p=",p,", q=",(p-bArg)/4);
		1/0;
	);
	if((p-bArg)%4==0,
		q=(p-bArg)/4;
		\\ abs(phi8)=1 means that (a,b)=(rho_{k-1}, rho_{k-1}-4*pi_{k})
		if(abs(phi8)==1,
			nArray=get_rho_index(a);
			if(dbg!=0,
				print("in check_a_b(): a=",a,", aArg=",aArg,", b=",b,", bArg=",bArg,", Phi_8=",phi8,", nArray=",nArray);
			);
			isFound=0;
			
			for(i=1,length(nArray),
				n=nArray[i];
				e=1;
				k=n+e;
				if(k<0 || k-e<0,
					printf("   WARNING: neg-index pi/rho not yet implemented: k=%2d, e=%2d: a=%9d, b=%9d, p=%9d, q=%9d, b=rho_n for n=%2d, Phi_8=%2d\n",k,e,a,b,p,q,n,phi8);
				);
				if(k>-1 && k-e>-1,
					b1a=get_rho(k-e)-4*get_pi(k);
					if(b1a==b,
						printf("FOUND rho: a=%9d, b=%9d, p=%9d, q=%9d, (a,b)=%s(rho_(k-e), rho_(k-e)-4pi_k ) for k=%2d, e=%2d, Phi_8=%2d\n",aArg,bArg,p,q,sgn,k,e,phi8);
						isFound=1;
					);
				);
				e=-1;
				k=n+e;
				if(k<0 || k-e<0,
					printf("   WARNING: neg-index pi/rho not yet implemented: k=%2d, e=%2d: a=%9d, b=%9d, p=%9d, q=%9d, a=2*pi_n for n=%2d, Phi_8=%2d\n",k,e,a,b,p,q,n,phi8);
				);
				if(k>-1 && k-e>-1,
					b1b=get_rho(k-e)-4*get_pi(k);
					if(b1b==b,
						printf("FOUND rho: a=%9d, b=%9d, p=%9d, q=%9d, (a,b)=%s(rho_(k-e), rho_(k-e)-4pi_k ) for k=%2d, e=%2d, Phi_8=%2d\n",aArg,bArg,p,q,sgn,k,e,phi8);
						isFound=1;
					);
				);
			);
			if(isFound==0,
				printf("\n(*) NOT FOUND b in rho: a=%9d, b=%9d, p=%9d, q=%9d, a= rho_n for n=%3d, b1a=%3d, b1b=%3d, Phi_8=%2d\n\n",a,b,p,q,n,b1a,b1b,phi8);
				1/0;
			);
			if(length(nArray)==0,
				printf("\n(*) NOT FOUND a in rho: a=%9d, b=%9d, p=%9d, q=%9d, Phi_8=%2d\n\n",a,b,p,q,phi8);
				1/0;
			);
		);
		
		\\ abs(phi8)=2 means it (a,b)=(2*pi_{k-1}, 2*pi_{k-1}-4*rho_{k})
		if(abs(phi8)==2,
			n=get_pi_index(a/2);
			if(n!=-1000,
				isFound=0;
				e=1;
				k=n+e;
				if(k<0 || k-e<0,
					printf("   WARNING: neg-index pi not yet implemented: k=%2d, e=%2d: a=%9d, b=%9d, p=%9d, q=%9d, a=2*pi_n for n=%2d, Phi_8=%2d\n",k,e,a,b,p,q,n,phi8);
				);
				if(k>-1 && k-e>-1,
					b1=2*get_pi(k-e)-4*get_rho(k);
					if(b1==b,
						printf("FOUND  pi: a=%9d, b=%9d, p=%9d, q=%9d, (a,b)=%s(2pi_(k-e), 2pi_(k-e)-4rho_k) for k=%2d, e=%2d, Phi_8=%2d\n",aArg,bArg,p,q,sgn,k,e,phi8);
						isFound=1;
					);
				);
				e=-1;
				k=n+e;
				if(k<0 || k-e<0,
					printf("   WARNING: neg-index pi not yet implemented: k=%2d, e=%2d: a=%9d, b=%9d, p=%9d, q=%9d, a=2*pi_n for n=%2d, Phi_8=%2d\n",k,e,a,b,p,q,n,phi8);
				);
				if(k>-1 && k-e>-1,
					b1=2*get_pi(k-e)-4*get_rho(k);
					if(b1==b,
						printf("FOUND  pi: a=%9d, b=%9d, p=%9d, q=%9d, (a,b)=%s(2pi_(k-e), 2pi_(k-e)-4rho_k) for k=%2d, e=%2d, Phi_8=%2d\n",aArg,bArg,p,q,sgn,k,e,phi8);
						isFound=1;
					);
				);
				if(isFound==0,
					printf("\n(*) NOT FOUND b in pi: a=%9d, b=%9d, p=%9d, q=%9d, a=2*pi_n for n=%3d, b1a=%2d, b1b=%2d, Phi_8=%2d\n",a,b,p,q,n,b1a,b1b,phi8);
					1/0;
				);
			);
			if(n==-1000,
				printf("\n(*) NOT FOUND a/2 in pi: a=%9d, b=%9d, p=%9d, q=%9d, Phi_8=%2d\n",a,b,p,q,phi8);
				1/0;
			);
		);
	);
}

\\ this is the sequence rho in Abouzaid's paper
\\ returns [] if not found
\\ only looks for non-negative indices
\\ 7 May 2024
get_rho_index(f)={
	my(i);
	
	if(f==1,return([0,1]));
	i=get_index(f,1,1,2,1,"rho");
	if(i!=-1000,
		return([i]);
	);
	return([]);
}

\\ this is the sequence pi in Abouzaid's paper
\\ returns -1000 if not found
\\ only looks for non-negative indices
\\ 7 May 2024
get_pi_index(f)={
	return(get_index(f,0,1,2,1,"pi"));
}

\\ this is the sequence rho in Abouzaid's paper
\\ 7 May 2024
get_rho(k)={
	return(get_kth_element(k,1,1,2,1,"rho"));
}

\\ just for testing/checking the values of the sequence
\\ 7 May 2024
output_rho()={
	output_sequence(0,1,2,1);
}

\\ this is the sequence pi in Abouzaid's paper
\\ 7 May 2024
get_pi(k)={
	return(get_kth_element(k,0,1,2,1,"pi"));
}

\\ just for testing/checking the values of the sequence
\\ 7 May 2024
output_pi()={
	output_sequence(0,1,2,1);
}
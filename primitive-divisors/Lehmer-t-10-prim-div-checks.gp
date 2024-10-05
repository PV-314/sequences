\\ \r sequences\primitive-divisors\Lehmer-t-10-prim-div-checks.gp

read("sequences\\primitive-divisors\\Lehmer-utils.gp");

\\ use \p 50 (or higher)
\\ this, and tol=10.0^(-30), will work with bnd up to 10^8

\\ entry point. Call this function. From outside, the other functions can be ignored.
t10_check(dbg=0)={
	my(bnd);
	
	bnd=10^7;
	for (a=-bnd,bnd,
		check_a_c(a,1,dbg);
		check_a_c(a,-1,dbg);
		check_a_c(a,5,dbg);
		check_a_c(a,-5,dbg);
	);
}

\\ the expressions for b come from f:=a*a+10*a*b+5*b*b:solve(f=16*c,b);
\\ c means the constant on the right-hand side. It will be \pm 1, \pm 5 here
check_a_c(a,c,dbg=0)={
	my(b);
	if(issquare(5*a*a+20*c),
		b=-a+2*sqrtint(5*a*a+20*c)/5;
		if(denominator(b)==1,
			check_alpha_beta(a,b,dbg);
			if(5*a*a+20*c!=0,
				b=-a-2*sqrtint(5*a*a+20*c)/5;
				check_alpha_beta(a,b,dbg);
			);
		);
	);
}

\\ make alpha and beta from a and b
check_alpha_beta(a,b,dbg=0)={
	my(al,be,phi10,tol);
	
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
		phi10=(a*a+10*a*b+5*b*b)/16;
		if(dbg!=0,
			print("in check_alpha_beta(): checking Phi_10=",phi10," for a=",a,", b=",b);
		);
		if(abs(phi10)==1 || abs(phi10)==5,
			if(dbg!=0,
				print("in check_alpha_beta(): checking sequence for a=",a,", b=",b,", Phi_10=",phi10);
			);
			check_a_b(a,b,phi10,dbg);
		);
	);
}

\\ do a and b come from where I want/expect?
check_a_b(aArg,bArg,phi10,dbg=0)={
	my(a,a1a,a1b,b,e,isFound,k,kMin,n,nArray,p,q,sgn);
	
	kMin=-2;
	a=aArg;
	b=bArg;
	sgn="+";
	if(bArg<0,
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
		if(abs(phi10)==1,
			nArray=get_fibonacci_index(b);
			if(dbg!=0,
				print("in check_a_b(): a=",a,", aArg=",aArg,", b=",b,", bArg=",bArg,", Phi_10=",phi10,", nArray=",nArray);
			);
			isFound=0;
			for(i=1,length(nArray),
				n=nArray[i];
				e=1;
				k=n+2*e;
				if(k<kMin || k-2*e<kMin,
					printf("   WARNING: index<%2d Fibon not yet implemented: k=%2d, e=%2d, k-2e=%2d: a=%8d, b=%10d, p=%8d, q=%9d, a=F_n for n=%2d, Phi_10=%2d\n",kMin,k,e,k-2*e,a,b,p,q,n,phi10);
				);
				if(k>kMin-1 && k-2*e>kMin-1,
					a1a=fibonacci(k-2*e)-4*fibonacci(k);
					if(a1a==a,
						printf("FOUND Fibon: a=%8d, b=%10d, p=%8d, q=%9d, (a,b)=%s(phi_(k-2e)-4phi_k, phi_(k-2e)) for k=%2d, e=%2d, Phi_10=%2d\n",a,b,p,q,sgn,k,e,phi10);
						isFound=1;
					);
				);
				e=-1;
				k=n+2*e;
				if(k<kMin || k-2*e<kMin,
					printf("   WARNING: index<%2d Fibon not yet implemented: k=%2d, e=%2d, k-2e=%2d: a=%8d, b=%10d, p=%8d, q=%9d, a=F_n for n=%2d, Phi_10=%2d\n",kMin,k,e,k-2*e,a,b,p,q,n,phi10);
				);
				if(k>kMin-1 && k-2*e>kMin-1,
					a1b=fibonacci(k-2*e)-4*fibonacci(k);
					if(a1b==a,
						printf("FOUND Fibon: a=%8d, b=%10d, p=%8d, q=%9d, (a,b)=%s(phi_(k-2e)-4phi_k, phi_(k-2e)) for k=%2d, e=%2d, Phi_10=%2d\n",a,b,p,q,sgn,k,e,phi10);
						isFound=1;
					);
				);
			);
			if(isFound==0,
				printf("\n(*) NOT FOUND a in Fibon: a=%8d, b=%10d, p=%8d, q=%9d, a=F_n for n=%3d, a1a=%3d, a1b=%3d, Phi_10=%2d\n\n",a,b,p,q,n,a1a,a1b,phi10);
				1/0;
			);
			if(length(nArray)==0,
				printf("\n(*) NOT FOUND b in Fibon: a=%8d, b=%10d, p=%8d, q=%9d, Phi_10=%2d\n\n",a,b,p,q,phi10);
				1/0;
			);
		);
		if(abs(phi10)==5,
			nArray=get_lucas_index(b);
			isFound=0;
			for(i=1,length(nArray),
				n=nArray[i];
				e=1;
				k=n+2*e;
				if(k<kMin || k-2*e<kMin,
					printf("   WARNING: index<%2d Lucas not yet implemented: k=%2d, e=%2d: k-2e=%2d, a=%8d, b=%10d, p=%8d, q=%9d, a=L_n for n=%2d, Phi_10=%2d\n",kMin,k,e,k-2*e,a,b,p,q,n,phi10);
				);
				if(k>kMin-1 && k-2*e>kMin-1,
					a1a=get_lucas(k-2*e)-4*get_lucas(k);
					if(a1a==a,
						printf("FOUND Lucas: a=%8d, b=%10d, p=%8d, q=%9d, (a,b)=%s(psi_(k-2e)-4psi_k, psi_(k-2e)) for k=%2d, e=%2d, Phi_10=%2d\n",aArg,bArg,p,q,sgn,k,e,phi10);
						isFound=1;
					);
				);
				e=-1;
				k=n+2*e;
				if(k<kMin || k-2*e<kMin,
					printf("   WARNING: index<%2d Lucas not yet implemented: k=%2d, e=%2d, k-2e=%2d: a=%8d, b=%10d, p=%8d, q=%9d, a=L_n for n=%2d, Phi_10=%2d\n",kMin,k,e,k-2*e,a,b,p,q,n,phi10);
				);
				if(k>kMin-1 && k-2*e>kMin-1,
					a1b=get_lucas(k-2*e)-4*get_lucas(k);
					if(a1b==a,
						printf("FOUND Lucas: a=%8d, b=%10d, p=%8d, q=%9d, (a,b)=%s(psi_(k-2e)-4psi_k, psi_(k-2e)) for k=%2d, e=%2d, Phi_10=%2d\n",aArg,bArg,p,q,sgn,k,e,phi10);
						isFound=1;
					);
				);
			);
			if(isFound==0,
				printf("\n(*) NOT FOUND a in Lucas: a=%8d, b=%10d, p=%8d, q=%9d, a=L_n for n=%3d, a1a=%2d, a1b=%2d, k=%2d, e=%2d, Phi_10=%2d\n\n",a,b,p,q,n,a1a,a1b,k,e,phi10);
				1/0;
			);
			if(length(nArray)==0,
				printf("\n(*) NOT FOUND b in Lucas: a=%8d, b=%10d, p=%8d, q=%9d, Phi_10=%2d\n\n",a,b,p,q,phi10);
				1/0;
			);
		);
	);
}

\\ returns [] if not found
\\ only looks for indices \geq -2
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

output_lucas_sequence()={
	output_sequence(2,1,1,1);
}

get_lucas(k)={
	if(k==-2,return(3));
	if(k==-1,return(-1));
	return(get_kth_element(k,2,1,1,1,"Lucas(psi)"));
}

\\ returns [] if not found
\\ only looks for indices \geq -2
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
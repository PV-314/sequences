\\ \r sequences\primitive-divisors\Lehmer-t-12-prim-div-checks.gp

read("sequences\\primitive-divisors\\Lehmer-utils.gp");

\\ uses Abouzaid indexing

\\ use \p 50 (or higher)
\\ this, and tol=10.0^(-30), will work with bnd up to 10^8

\\ entry point. Call this function. From outside, the other functions can be ignored.
t12_check(dbg=0)={
	my(bnd);
	
	bnd=10^7;
	for (a=-bnd,bnd,
		check_a_c(a,1,dbg);
		check_a_c(a,-1,dbg);
		check_a_c(a,2,dbg);
		check_a_c(a,-2,dbg);
		check_a_c(a,3,dbg);
		check_a_c(a,-3,dbg);
		check_a_c(a,6,dbg);
		check_a_c(a,-6,dbg);
	);
}

\\ the expressions for b come from phi12:=a*a+14*a*b+b*b:solve(phi12=16*c,b);
check_a_c(a,c,dbg=0)={
	my(b);
	
	if(issquare(3*a*a+c),
		b=-7*a+4*sqrtint(3*a^2+c);
		check_alpha_beta(a,b,dbg);
		if(3*a*a+c!=0,
			b=-7*a-4*sqrtint(3*a^2+c);
			check_alpha_beta(a,b,dbg);
		);
	);
}

check_alpha_beta(a,b,dbg=0)={
	my(al,be,phi12,tol);
	
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
		phi12=(a*a+14*a*b+b*b)/16;
		if(dbg!=0,
			print("in check_alpha_beta(): checking Phi_12=",phi12," for a=",a,", b=",b);
		);
		if(abs(phi12)==1 || abs(phi12)==2 || abs(phi12)==3 || abs(phi12)==6,
			if(dbg!=0,
				print("in check_alpha_beta(): checking sequence for a=",a,", b=",b,", Phi_12=",phi12);
			);
			check_a_b(a,b,phi12,dbg);
		);
	);
}

check_a_b(a,b,phi12,dbg=0)={
	my(p,q);
	
	p=a;
	if((p-b)%4!=0,
		print("BAD: in in check_a_b(), a=",a,", b=",b,", p=",p,", q=",(p-b)/4);
		1/0;
	);
	if((p-b)%4==0,
		q=(p-b)/4;
		check_Phi12(a,b,p,q,phi12,dbg);
	);
}

check_Phi12(aArg,bArg,p,q,phi12,dbg=0)={
	my(a,b,b1a,b1b,e,isFound,k,kMin,n,name,sgn);

	kMin=0;
	name=get_sequence_name(phi12);
	sgn="+";
	a=aArg;
	b=bArg;
	if(aArg<0,
		a=-aArg;
		b=-bArg;
		sgn="-";
	);
	
	if(abs(phi12)==6,
		\\ to handle k=-1
		kMin=-1;
		if(aArg==-1 && bArg==-5,
			a=aArg;
			b=bArg;
			sgn="+";
		);
		if(aArg==1 && bArg==5,
			a=-aArg;
			b=-bArg;
			sgn="-";
		);
	);
	
	n=get_sequence_index(a,phi12);
	if(n!=-1000,
		isFound=0;
		e=1;
		k=n+e;
		if(k+e<kMin,
			printf("   WARNING: neg-index %s not yet implemented: k=%2d, e=%2d: a=%6d, b=%6d, p=%6d, q=%6d, a=zeta0_n for n=%2d, Phi_12=%2d\n",name,k,e,a,b,p,q,n,phi12);
		);
		if(k+e>kMin-1,
			b1a=get_sequence_value(k+e,phi12);
			if(b1a==-b,
				printf("FOUND %s: a=%8d, b=%10d, p=%8d, q=%9d, (a,b)=%s(%s_(k-e), -%s_(k+e)) for k=%2d, e=%2d, Phi_12=%2d\n",name,aArg,bArg,p,q,sgn,name,name,k,e,phi12);
				isFound=1;
			);
		);
		e=-1;
		k=n+e;
		if(k+e<kMin,
			printf("   WARNING: neg-index %s not yet implemented: k=%2d, e=%2d: a=%6d, b=%6d, p=%6d, q=%6d, b=%s_k for k=%2d, Phi_12=%2d\n",name,k,e,aArg,bArg,p,q,name,k,phi12);
		);
		if(k+e>kMin-1,
			b1b=get_sequence_value(k+e,phi12);
			if(b1b==-b,
				printf("FOUND %s: a=%8d, b=%10d, p=%8d, q=%9d, (a,b)=%s(%s_(k-e), -%s_(k+e)) for k=%2d, e=%2d, Phi_12=%2d\n",name,aArg,bArg,p,q,sgn,name,name,k,e,phi12);
				isFound=1;
			);
		);
		if(dbg!=0,
			printf("%s: aArg=%8d, a=%8d, bArg=%10d, b=%10d, p=%8d, q=%9d, k=%2d, e=%2d, Phi_12=%2d\n",name,aArg,a,bArg,b,p,q,k,e,phi12);
		);
		if(isFound==0,
			printf("NOT FOUND -b in %s: aArg=%6d, a=%6d, bArg=%6d, b=%6d, p=%6d, q=%6d, b=%s_k for k=%3d, b1a=%2d, b1b=%2d, Phi_12=%2d\n",name,aArg,a,bArg,b,p,q,name,k,b1a,b1b,phi12);
			1/0;
		);
	);
	if(n==-1000,
		printf("NOT FOUND  a in %s: a=%6d, b=%6d, p=%6d, q=%6d, Phi_12=%2d\n",name,aArg,bArg,p,q,phi12);
		1/0;
	);
}

get_sequence_index(a,phi12)={
	if(abs(phi12)==1,
		return(get_zeta0_index(a));
	);
	if(abs(phi12)==2,
		return(get_zeta2P1_index(a));
	);
	if(abs(phi12)==3,
		return(get_zeta1_index(a));
	);
	if(abs(phi12)==6,
		return(get_zeta3P1_index(a));
	);
	print("ERROR, bad Phi_12 in get_sequence_index() with Phi_12=",phi12);
	return(1/0);
}

get_sequence_name(phi12)={
	if(abs(phi12)==1,
		return("  zeta0");
	);
	if(abs(phi12)==2,
		return("zeta2P1");
	);
	if(abs(phi12)==3,
		return("  zeta1");
	);
	if(abs(phi12)==6,
		return("zeta3P1");
	);
	print("ERROR, bad Phi_12 in get_sequence_name() with Phi_12=",phi12);
	return(1/0);
}

get_sequence_value(k,phi12)={
	if(abs(phi12)==1,
		return(get_zeta0(k));
	);
	if(abs(phi12)==2,
		return(get_zeta2P1(k));
	);
	if(abs(phi12)==3,
		return(get_zeta1(k));
	);
	if(abs(phi12)==6,
		return(get_zeta3P1(k));
	);
	print("ERROR, bad Phi_12 in get_sequence_value() with Phi_12=",phi12);
	return(1/0);
}

\\ returns -1000 if not found
\\ only looks for non-negative indices
get_zeta0_index(f)={
	return(get_index(f,0,1,4,-1,"zeta0"));
}

get_zeta0(k)={
	return(get_kth_element(k,0,1,4,-1,"zeta0"));
}

\\ just for testing/checking the values of the sequence
output_zeta0()={
	output_sequence(0,1,4,-1);
}

get_zeta1(k)={
	return(get_kth_element(k,1,2,4,-1,"zeta1"));
}

\\ returns -1000 if not found
\\ only looks for non-negative indices
get_zeta1_index(f)={
	return(get_index(f,1,2,4,-1,"zeta1"));
}

\\ just for testing/checking the values of the sequence
output_zeta1()={
	output_sequence(1,2,4,-1);
}

get_zeta2P1(k)={
	my(e);
	e=1;
	return(get_kth_element(k,e,2*e+1,4,-1,"zeta2P1"));
}

\\ returns -1000 if not found
\\ only looks for non-negative indices
get_zeta2P1_index(f)={
	my(e);
	e=1;
	return(get_index(f,e,2*e+1,4,-1,"zeta2P1"));
}

\\ just for testing/checking the values of the sequence
output_zeta2P1()={
	my(e);
	e=1;
	output_sequence(e,2*e+1,4,-1);
}

get_zeta3P1(k)={
	my(e);

	\\ my extension for k \geq -1
	if(k==-1,return(-1));
	e=1;
	return(get_kth_element(k,1,3*e+2,4,-1,"zeta3P1"));
}

\\ returns -1000 if not found
\\ only looks for indices \geq -1
get_zeta3P1_index(f)={
	my(e);
	
	\\ my extension for k \geq -1
	if(f==-1,return(-1));
	e=1;
	return(get_index(f,1,3*e+2,4,-1,"zeta3P1"));
}

\\ just for testing/checking the values of the sequence
output_zeta3P1()={
	my(e);
	e=1;
	output_sequence(1,3*e+2,4,-1);
}
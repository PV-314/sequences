get_un(a,b,n,i)={
	my(al,be,tol,un);
	
	tol=10.0^(-8);
	al=(sqrt(a)+sqrt(b))/2;
	be=(sqrt(a)-sqrt(b))/2;
	un=(al^n-be^n)/(al-be);
	if(n%2==0,
		un=un/(al+be);
	);
	if(abs(imag(un))>tol,
		print("ERROR: NOT REAL u_",n,"=",un," for i=",i,", a=",a,", b=",b);
		1/0;
	);
	un=real(un);
	if(abs(un-round(un))>tol,
		print("ERROR: NOT INTEGER u_",n,"=",un," for i=",i,", a=",a,", b=",b);
		1/0;
	);
	un=round(un);
	return(un);
}

output_sequence(u0,u1,a1,a0,name="")={
	for(k=0,20,
		if(name=="",
			printf("n=%3d, u_n=%6d\n",k,get_kth_element(k,u0,u1,a1,a0,""));
		);
		if(name!="",
			printf("n=%3d, %s_n=%6d\n",k,name,get_kth_element(k,u0,u1,a1,a0,""));
		);
	);
}

\\ using recurrence u_{n+1}=a1*u_{n} + a0*u_{n-1} where u_0=u0 and u_1=u1
get_kth_element(k,u0,u1,a1,a0,name)={
	my(fC,fP,fT);
	
	if(k<0,
		printf("   %s for k=%3d: not implemented\n",name,k);
		return(-100000000000000);
	);
	if(k==0,return(u0));
	if(k==1,return(u1));
	fP=u0;
	fC=u1;
	for(n=2,k,
		fT=a1*fC+a0*fP;
		fP=fC;
		fC=fT;
	);
	return(fC);
}

\\ returns -1000 if not found
\\ only looks for non-negative indices
get_index(f,u0,u1,a1,a0,name="")={
	my(fC,fP,fT);
	
	if(f==u0,return(0));
	if(f==u1,return(1));
	fP=u0;
	fC=u1;
	for(k=2,200,
		fT=a1*fC+a0*fP;
		fP=fC;
		fC=fT;
		if(f==fC,
			return(k);
		);
	);
	printf("   %4d not found in %s\n",f,name);
	return(-1000);
}
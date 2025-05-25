\\ \r sequences\primitive-divisors\Lucas-t-2-sequences.gp

t2_check()={
	my(a,aBnd,b,d,isDegen,m,p,q,rtRat,sqrtD,tol,twoP,u2,u2P);
	
	tol=0.0000001;
	twoP=65536*65536;
	\\ the companion polynomial is x^2-a1*x+a0
	aBnd=3000;
	for(a1=-aBnd,aBnd,
	for(a0=-aBnd,aBnd,
		if(a0*a1!=0 && gcd(a0,a1)==1,
			q=a0;
			p=a1*a1;
			m=a1;
			if(m>0 && gcd(m,q)==1,
				d=a1*a1-4*a0;
				sqrtD=sqrt(d);
				rtRat=(a1+sqrtD)/(a1-sqrtD);
				u2=a1;
				a=m;
				b=p-4*q;
				isDegen=0;
				if(abs(rtRat^24-1)<tol,
					printf("degenerate: a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u2=%4d, rtRat=%9.6f+(%9.6f)*I: abs(rtRat^24-1)=%9.6e\n",a1,a0,a,b,m,p,q,u2,real(rtRat),imag(rtRat),abs(rtRat^24-1));
					isDegen=1;
				);
				if(isDegen==0,
					u2P=u2;
					while(gcd(u2P,d)!=1,
						u2P=u2P/gcd(u2P,d);
					);
					if(m==1 && abs(u2P)>1,
						printf("BAD, should be defective: a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u2=%4d, u2P=%4d\n",a1,a0,a,b,m,p,q,u2,u2P);
					);
					if(gcd(m,twoP)==m && q%2==1 && abs(u2P)>1,
						printf("BAD, should be defective: a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u2=%4d, u2P=%4d\n",a1,a0,a,b,m,p,q,u2,u2P);
					);
					if(abs(u2P)<2 && !(m==1 || (gcd(m,twoP)==m && q%2==1)),
						printf("BAD, should not be defective: a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u2=%4d, u2P=%4d\n",a1,a0,a,b,m,p,q,u2,u2P);
					);
				);
			);
		);
	);
	);
}
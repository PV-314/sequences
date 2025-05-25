\\ \r sequences\primitive-divisors\Lucas-t-3-sequences.gp

t3_check()={
	my(a,aBnd,b,b1,b2,d,isDegen,m,p,q,rtRat,sqrtD,threeP,tol,twoP,u1,u2,u3,u3P);
	
	tol=0.0000001;
	twoP=65536*65536;
	threeP=729*729*729*729;
	\\ the companion polynomial is x^2-a1*x+a0
	aBnd=3000;
	for(a1=-aBnd,aBnd,
	for(a0=-aBnd,aBnd,
		if(a0*a1!=0 && gcd(a0,a1)==1,
			m=a1;
			q=a0;
			if(m>0 && gcd(m,q)==1,
				d=a1*a1-4*a0;
				sqrtD=sqrt(d);
				rtRat=(a1+sqrtD)/(a1-sqrtD);
				p=a1*a1;
				a=m;
				b=p-4*q;
				u1=1;
				u2=a1;
				u3=a1*a1-a0;
				isDegen=0;
				if(abs(rtRat^24-1)<tol,
					printf("degenerate: a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u3=%4d, rtRat=%9.6f+(%9.6f)*I, abs(rtRat^24-1)=%9.6e\n",a1,a0,a,b,m,p,q,u3,real(rtRat),imag(rtRat),abs(rtRat^24-1));
					isDegen=1;
				);
				if(isDegen==0,
					u3P=u3;
					while(gcd(u3P,d*u2)!=1,
						u3P=u3P/gcd(u3P,d*u2);
					);
					if(abs(u3P)<2,
						isOK1=0;
						if(m>1 && b==4*1-3*m*m,
							isOK1=1;
						);
						if(b==4*(-1)-3*m*m,
							isOK1=1;
						);
						b1=(b+3*m*m)/4;
						b2=(b+3*m*m)/(-4);
						\\ need next check before k is a positive integer in this case in BHV:
						if(b1%3==0 || b2%3==0,
							b1=b1/gcd(b1,threeP);
							b2=b2/gcd(b2,threeP);
							if(m%3!=0 && (abs(b1)==1 || abs(b2)==1),
								isOK1=1;
							);
						);
						if(isOK1==0,
							printf("BAD, defective but not included: a1=%6d, a0=%6d, a=%6d, b=%6d, m=%6d, p=%6d, q=%6d, u3=%6d, u3P=%6d\n",a1,a0,a,b,m,p,q,u3,u3P);
						);
					);
					if(m>1 && (b==4*1-3*m*m || b==4*(-1)-3*m*m),
						if(abs(u3P)>1,
							printf("BAD, included but not defective(1): a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u3=%4d, u3P=%4d\n",a1,a0,a,b,m,p,q,u3,u3P);
						);
					);
					b1=(b+3*m*m)/4;
					b2=(b+3*m*m)/(-4);
					if(m%3!=0 && (b1%3==0 || b2%3==0),
						b1=b1/gcd(b1,threeP);
						b2=b2/gcd(b2,threeP);
						if(abs(b1)==1 || abs(b2)==1,
							if(abs(u3P)>1,
								printf("BAD, included but not defective(2): a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u3=%4d, u3P=%4d\n",a1,a0,a,b,m,p,q,u3,u3P);
							);
						);
					);
				);
			);
		);
	);
	);
}
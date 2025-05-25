\\ \r sequences\primitive-divisors\Lucas-t-4-sequences.gp

t4_check()={
	my(a,aBnd,b,b1,b2,d,isDegen,m,p,q,rtRat,sqrtD,threeP,tol,twoP,u1,u2,u3,u4,u4P);
	
	tol=0.0000001;
	twoP=65536*65536;
	threeP=729*729*729*729;
	\\ the companion polynomial is x^2-a1*x+a0
	aBnd=3000;
	for(a1=-aBnd,aBnd,
	for(a0=-aBnd,aBnd,
			m=a1;
			q=a0;
			p=a1*a1;
			a=m;
			b=p-4*q;
			if(a==1 && b==1,
				printf("(1): a1=%6d, a0=%6d, a=%6d, b=%6d, m=%6d, p=%6d, q=%6d\n",a1,a0,a,b,m,p,q);
			);
		if(a0*a1!=0 && gcd(a0,a1)==1,
			m=a1;
			q=a0;
			p=a1*a1;
			a=m;
			b=p-4*q;
			if(a==1 && b==1,
				printf("(2): a1=%6d, a0=%6d, a=%6d, b=%6d, m=%6d, p=%6d, q=%6d\n",a1,a0,a,b,m,p,q);
			);
			if(m>0 && gcd(m,q)==1,
				d=a1*a1-4*a0;
				sqrtD=sqrt(d);
				rtRat=(a1+sqrtD)/(a1-sqrtD);
				u1=1;
				u2=a1;
				u3=a1*a1-a0;
				u4=a1*u3-a0*u2;
				isDegen=0;
				if(abs(rtRat^24-1)<tol,
					printf("degenerate: a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u4=%4d, rtRat=%9.6f+(%9.6f)*I, abs(rtRat^24-1)=%9.6e\n",a1,a0,a,b,m,p,q,u4,real(rtRat),imag(rtRat),abs(rtRat^24-1));
					isDegen=1;
				);
				if(isDegen==0,
					if(a==1 && b==1,
						printf("(3): a1=%6d, a0=%6d, a=%6d, b=%6d, m=%6d, p=%6d, q=%6d, u2=%6d, u3=%6d, u4=%6d\n",a1,a0,a,b,m,p,q,u2,u3,u4);
					);
					u4P=u4;
					while(gcd(u4P,d*u2*u3)!=1,
						u4P=u4P/gcd(u4P,d*u2*u3);
					);
					if(abs(u4P)<2,
						isOK1=0;
						if(m>1 && m%2==1 && (b==2-m*m || b==-2-m*m),
							isOK1=1;
						);
						if(m>1 && m%2==0 && b==-4-m*m,
							isOK1=1;
						);
						if(m>2 && m%2==0 && b==4-m*m,
							isOK1=1;
						);
						if(isOK1==0,
							printf("BAD, defective but not included: a1=%6d, a0=%6d, a=%6d, b=%6d, m=%6d, p=%6d, q=%6d, u4=%6d, u4P=%6d\n",a1,a0,a,b,m,p,q,u4,u4P);
						);
					);
					if(m>1 && m%2==1 && (b==2-m*m || b==-2-m*m),
						if(abs(u4P)>1,
							printf("BAD, included but not defective(1): a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u4=%4d, u4P=%4d\n",a1,a0,a,b,m,p,q,u4,u4P);
						);
					);
					if(m>1 && m%2==0 && b==-4-m*m,
						if(abs(u4P)>1,
							printf("BAD, included but not defective(1): a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u4=%4d, u4P=%4d\n",a1,a0,a,b,m,p,q,u4,u4P);
						);
					);
					if(m>2 && m%2==0 && b==4-m*m,
						if(abs(u4P)>1,
							printf("BAD, included but not defective(1): a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u4=%4d, u4P=%4d\n",a1,a0,a,b,m,p,q,u4,u4P);
						);
					);
				);
			);
		);
	);
	);
}
\\ \r sequences\primitive-divisors\Lucas-t-6-sequences.gp

t6_check()={
	my(a,aBnd,b,b1,b2,d,isDegen,m,p,q,rtRat,sqrtD,threeP,tol,twoP,u1,u2,u3,u4,u5,u6,u6P);
	
	tol=0.0000001;
	twoP=65536*65536*65536;
	threeP=729*729*729*729;
	\\ the companion polynomial is x^2-a1*x+a0
	aBnd=3000;
	for(a1=-aBnd,aBnd,
	for(a0=-aBnd,aBnd,
		if(a0*a1!=0 && gcd(a0,a1)==1,
			m=a1;
			q=a0;
			p=a1*a1;
			a=m;
			b=p-4*q;
			if(m>0 && gcd(m,q)==1,
				d=a1*a1-4*a0;
				sqrtD=sqrt(d);
				rtRat=(a1+sqrtD)/(a1-sqrtD);
				u1=1;
				u2=a1;
				u3=a1*a1-a0;
				u4=a1*u3-a0*u2;
				u5=a1*u4-a0*u3;
				u6=a1*u5-a0*u4;
				isDegen=0;
				if(abs(rtRat^24-1)<tol,
					printf("degenerate: a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u4=%4d, rtRat=%9.6f+(%9.6f)*I, abs(rtRat^24-1)=%9.6e\n",a1,a0,a,b,m,p,q,u4,real(rtRat),imag(rtRat),abs(rtRat^24-1));
					isDegen=1;
				);
				if(isDegen==0,
					u6P=u6;
					while(gcd(u6P,d*u2*u3*u4*u5)!=1,
						u6P=u6P/gcd(u6P,d*u2*u3*u4*u5);
					);
					if(abs(u6P)<2,
						isOK1=0;
						if(m>3 && m%3!=0 && b==(4-m*m)/3,
							isOK1=1;
						);
						if(m%3==0 && b==4-m*m/3,
							isOK1=1;
						);
						if(m%3==0 && b==-4-m*m/3,
							isOK1=1;
						);
						b1=3*b+m*m;
						if((m%6==1 || m%6==5) && (b1%8==0 && abs(b1/gcd(twoP,b1))==1),
							eB=log(abs(b1))/log(2);
							if(abs(eB-round(eB))>0.01,
								printf("BAD eB: a1=%6d, a0=%6d, a=%6d, b=%6d, m=%6d, p=%6d, q=%6d, b1=%6d, eB=%9.6f, u6=%6d, u6P=%6d\n",a1,a0,a,b,m,p,q,b1,eB,u6,u6P);
							);
							eB=round(eB);
							if(eB%2==0 && b1>0,
								isOK1=1;
							);
							if(eB%2!=0 && b1<0,
								isOK1=1;
							);
						);
						b2=b+m*m/3;
						if(m%6==3 && (b2%8==0 && abs(b2/gcd(twoP,b2))==1),
							isOK1=1;
						);
						if(isOK1==0,
							printf("BAD, defective but not included: a1=%6d, a0=%6d, a=%6d, b=%6d, m=%6d, p=%6d, q=%6d, u6=%6d, u6P=%6d\n",a1,a0,a,b,m,p,q,u6,u6P);
						);
					);
					if(m>3 && m%3!=0 && b==(4-m*m)/3,
						if(abs(u6P)>1,
							printf("BAD, included but not defective(1): a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u6=%4d, u6P=%4d\n",a1,a0,a,b,m,p,q,u6,u6P);
						);
					);
					if(m%3==0 && b==4-m*m/3,
						if(abs(u6P)>1,
							printf("BAD, included but not defective(2): a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u6=%4d, u6P=%4d\n",a1,a0,a,b,m,p,q,u6,u6P);
						);
					);
					if(m%3==0 && b==-4-m*m/3,
						if(abs(u6P)>1,
							printf("BAD, included but not defective(3): a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u6=%4d, u6P=%4d\n",a1,a0,a,b,m,p,q,u6,u6P);
						);
					);
					b1=3*b+m*m;
					if((m%6==1 || m%6==5) && (b1%8==0 && abs(b1/gcd(twoP,b1))==1),
						eB=log(abs(b1))/log(2);
						if(abs(eB-round(eB))>0.01,
							printf("BAD eB: a1=%6d, a0=%6d, a=%6d, b=%6d, m=%6d, p=%6d, q=%6d, b1=%6d, eB=%9.6f, u6=%6d, u6P=%6d\n",a1,a0,a,b,m,p,q,b1,eB,u6,u6P);
						);
						eB=round(eB);
						if((eB%2==0 && b1>0) || (eB%2!=0 && b1<0),
							if(abs(u6P)>1,
								printf("BAD, included but not defective(4): a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u6=%4d, u6P=%4d\n",a1,a0,a,b,m,p,q,u6,u6P);
							);
						);
					);
					b2=b+m*m/3;
					if(m%6==3 && (b2%8==0 && abs(b2/gcd(twoP,b2))==1),
						if(abs(u6P)>1,
							printf("BAD, included but not defective(5): a1=%4d, a0=%4d, a=%4d, b=%4d, m=%4d, p=%4d, q=%4d, u6=%4d, u6P=%4d\n",a1,a0,a,b,m,p,q,u6,u6P);
						);
					);
				);
			);
		);
	);
	);
}
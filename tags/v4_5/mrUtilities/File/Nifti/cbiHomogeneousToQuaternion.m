function [quatern, qfac, pixdim]=cbiHomogeneousToQuaternion( R )
% Calculates quaternion parameters from a 4x4 matrix
% [quatern, qfac, pixdim]=cbiHomogeneousToQuaternion( M44 )
% 
% - OR - (using a 3x3 rotation matrix)
% [quatern, qfac, pixdim]=cbiHomogeneousToQuaternion( R )
% 

% check args
    
  if (nargin~=1)
    error('Wrong number of input arguments!');
  end

  if (size(R)==[4,4] | size(R)==[3,3])
    R=R(1:3,1:3);
  else
    error('Input argument must be a homogeneous 4x4 matrix or 3x3 rotation matrix');
  end
  
  if (nargout~=3)
    error('Wrong number of output arguments!');
  end

  [quatern,qfac, pixdim]=calcQuatFromRotMat(R);

  return
  
function [quatern,qfac,pixdim]=calcQuatFromRotMat(R)
  % BEGIN calc
  % assume R set
  
  %  /* compute lengths of each column; these determine grid spacings  */
 

  xd = sqrt( sum(R(:,1).^2) );
  yd = sqrt( sum(R(:,2).^2) );
  zd = sqrt( sum(R(:,3).^2) );
 
  
% $$$    /* if a column length is zero, patch the trouble */
  
  if ( xd == 0 )
    R(1,1) = 1 ; R(2,1) = 0; R(3,1) = 0 ; xd = 1 ;
  end
  if ( yd == 0 )
    R(2,2) = 1 ; R(1,2) = 0; R(3,2) = 0 ; yd = 1 ; 
  end
  if ( zd == 0 )
    R(3,3) = 1 ; R(1,3) = 0; R(2,3) = 0 ; zd = 1 ; 
  end
  pixdim=[xd yd zd];
% $$$    /* normalize the columns */
    
  R(:,1)=R(:,1)./xd;
  R(:,2)=R(:,2)./yd;
  R(:,3)=R(:,3)./zd;
  
  R=polarDecomp33(R);
  zd=det(R);
  
  if (zd>0)
    qfac=1;
  else
    qfac=-1;
    R(:,3)=-R(:,3);
  end
  
  a=R(1,1)+R(2,2)+R(3,3)+1;
  if (a>0.5)
    a = 0.5 * sqrt(a) ;
    b = 0.25 * (R(3,2)-R(2,3)) / a ;
    c = 0.25 * (R(1,3)-R(3,1)) / a ;
    d = 0.25 * (R(2,1)-R(1,2)) / a ;
  else
    xd = 1.0 + R(1,1) - (R(2,2)+R(3,3)) ;
    yd = 1.0 + R(2,2) - (R(1,1)+R(3,3)) ;
    zd = 1.0 + R(3,3) - (R(1,1)+R(2,2)) ;
    if( xd > 1.0 )
      b = 0.5 * sqrt(xd) ;
      c = 0.25* (R(1,2)+R(2,1)) / b ;
      d = 0.25* (R(1,3)+R(3,1)) / b ;
      a = 0.25* (R(3,2)-R(2,3)) / b ;
    elseif ( yd > 1.0 )
      c = 0.5 * sqrt(yd) ;
      b = 0.25* (R(1,2)+R(2,1)) / c ;
      d = 0.25* (R(2,3)+R(3,2)) / c ;
      a = 0.25* (R(1,3)-R(3,1)) / c ;
    else
      d = 0.5 * sqrt(zd) ;
      b = 0.25* (R(1,3)+R(3,1)) / d ;
      c = 0.25* (R(2,3)+R(3,2)) / d ;
      a = 0.25* (R(2,1)-R(1,2)) / d ;
    end
    if( a < 0.0 )
      b=-b; c=-c; d=-d; a=-a; 
    end
    
  end
  quatern=[b c d];
  
  return;

  
  
function OM=polarDecomp33( M );
  
   X = M;

   gam = det(X) ;
   while ( gam == 0.0 )
     gam = 0.00001 * ( 0.001 + norm(X,'inf'));
     X(1,1)=X(1,1)+gam;
     X(2,2)=X(2,2)+gam;
     X(3,3)=X(3,3)+gam;
     gam = det(X);
   end
   dif=1;
   k=0;
   while (1)
     Y = inv(X);
     if( dif > 0.3 )
       alp = sqrt( norm(X,'inf') * norm(X','inf'));
       bet = sqrt( norm(Y,'inf') * norm(Y','inf'));
       gam = sqrt( bet / alp ) ;
       gmi = 1.0 / gam ;
     else
       gam = 1.0;
       gmi = 1.0 ;
     end

     Z(1,1) = 0.5 * ( gam*X(1,1) + gmi*Y(1,1) ) ;     
     Z(1,2) = 0.5 * ( gam*X(1,2) + gmi*Y(2,1) ) ;
     Z(1,3) = 0.5 * ( gam*X(1,3) + gmi*Y(3,1) ) ;
     Z(2,1) = 0.5 * ( gam*X(2,1) + gmi*Y(1,2) ) ;
     Z(2,2) = 0.5 * ( gam*X(2,2) + gmi*Y(2,2) ) ;
     Z(2,3) = 0.5 * ( gam*X(2,3) + gmi*Y(3,2) ) ;
     Z(3,1) = 0.5 * ( gam*X(3,1) + gmi*Y(1,3) ) ;
     Z(3,2) = 0.5 * ( gam*X(3,2) + gmi*Y(2,3) ) ;
     Z(3,3) = 0.5 * ( gam*X(3,3) + gmi*Y(3,3) ) ;

     dif = abs(Z(1,1)-X(1,1))+abs(Z(1,2)-X(1,2)) ...
	   +abs(Z(1,3)-X(1,3))+abs(Z(2,1)-X(2,1)) ...
	   +abs(Z(2,2)-X(2,2))+abs(Z(2,3)-X(2,3)) ...	   
	   +abs(Z(3,1)-X(3,1))+abs(Z(3,2)-X(3,2)) ...
	   +abs(Z(3,3)-X(3,3));

     k = k+1 ;
     if ( k > 100 | dif < 3.e-6 ) 
       break ;
     end
     X = Z ;
   end

   OM=Z;
   return


% 
% $$$ /* Given the 3x4 upper corner of the matrix R, compute the quaternion
% $$$        parameters that fit it.  See comments in nifti1.h for details.
% $$$      - Any NULL pointer on input won't get assigned (e.g., if you don't want
% $$$        dx,dy,dz, just pass NULL in for those pointers).
% $$$      - If the 3 input matrix columns are NOT orthogonal, they will be
% $$$        orthogonalized prior to calculating the parameters, using
% $$$        the polar decomposition to find the orthogonal matrix closest
% $$$        to the column-normalized input matrix.
% $$$      - However, if the 3 input matrix columns are NOT orthogonal, then
% $$$        the matrix produced by quatern_to_mat44 WILL have orthogonal
% $$$        columns, so it won't be the same as the matrix input here.
% $$$        This "feature" is because the NIFTI 'qform' transform is
% $$$        deliberately not fully general -- it is intended to model a volume
% $$$        with perpendicular axes.
% $$$      - If the 3 input matrix columns are not even linearly independent,
% $$$        you'll just have to take your luck, won't you?
% $$$ -----------------------------------------------------------------------------*/
% $$$ 
% $$$ void mat44_to_quatern( mat44 R ,
% $$$                        float *qb, float *qc, float *qd, QUATERN params
% $$$                        float *qx, float *qy, float *qz, qoffsets
% $$$                        float *dx, float *dy, float *dz, float *qfac ) pixdims
% $$$ {
% $$$    double r11,r12,r13 , r21,r22,r23 , r31,r32,r33 ;
% $$$    double xd,yd,zd , a,b,c,d ;
% $$$    mat33 P,Q ;
% $$$ 
% $$$    /* offset outputs are read write out of input matrix  */
% $$$ 
% $$$    ASSIF(qx,R.m[0][3]) ; ASSIF(qy,R.m[1][3]) ; ASSIF(qz,R.m[2][3]) ;
% $$$ 
% $$$    /* load 3x3 matrix into local variables */
% $$$ 
% $$$    r11 = R.m[0][0] ; r12 = R.m[0][1] ; r13 = R.m[0][2] ;
% $$$    r21 = R.m[1][0] ; r22 = R.m[1][1] ; r23 = R.m[1][2] ;
% $$$    r31 = R.m[2][0] ; r32 = R.m[2][1] ; r33 = R.m[2][2] ;
% $$$ 
% $$$    /* compute lengths of each column; these determine grid spacings  */
% $$$ 
% $$$    xd = sqrt( r11*r11 + r21*r21 + r31*r31 ) ;
% $$$    yd = sqrt( r12*r12 + r22*r22 + r32*r32 ) ;
% $$$    zd = sqrt( r13*r13 + r23*r23 + r33*r33 ) ;
% $$$ 
% $$$    /* if a column length is zero, patch the trouble */
% $$$ 
% $$$    if( xd == 0.0l ){ r11 = 1.0l ; r21 = r31 = 0.0l ; xd = 1.0l ; }
% $$$    if( yd == 0.0l ){ r22 = 1.0l ; r12 = r32 = 0.0l ; yd = 1.0l ; }
% $$$    if( zd == 0.0l ){ r33 = 1.0l ; r13 = r23 = 0.0l ; zd = 1.0l ; }
% $$$ 
% $$$    /* assign the output lengths */
% $$$ 
% $$$    ASSIF(dx,xd) ; ASSIF(dy,yd) ; ASSIF(dz,zd) ;
% $$$ 
% $$$    /* normalize the columns */
% $$$ 
% $$$    r11 /= xd ; r21 /= xd ; r31 /= xd ;
% $$$    r12 /= yd ; r22 /= yd ; r32 /= yd ;
% $$$    r13 /= zd ; r23 /= zd ; r33 /= zd ;
% $$$ 
% $$$    /* At this point, the matrix has normal columns, but we have to allow
% $$$       for the fact that the hideous user may not have given us a matrix
% $$$       with orthogonal columns.
% $$$ 
% $$$       So, now find the orthogonal matrix closest to the current matrix.
% $$$ 
% $$$       One reason for using the polar decomposition to get this
% $$$       orthogonal matrix, rather than just directly orthogonalizing
% $$$       the columns, is so that inputting the inverse matrix to R
% $$$       will result in the inverse orthogonal matrix at this point.
% $$$       If we just orthogonalized the columns, this wouldn't necessarily hold. */
% $$$ 
% $$$    Q.m[0][0] = r11 ; Q.m[0][1] = r12 ; Q.m[0][2] = r13 ; /* load Q */
% $$$    Q.m[1][0] = r21 ; Q.m[1][1] = r22 ; Q.m[1][2] = r23 ;
% $$$    Q.m[2][0] = r31 ; Q.m[2][1] = r32 ; Q.m[2][2] = r33 ;
% $$$ 
% $$$    P = mat33_polar(Q) ;  /* P is orthog matrix closest to Q */
% $$$ 
% $$$    r11 = P.m[0][0] ; r12 = P.m[0][1] ; r13 = P.m[0][2] ; /* unload */
% $$$    r21 = P.m[1][0] ; r22 = P.m[1][1] ; r23 = P.m[1][2] ;
% $$$    r31 = P.m[2][0] ; r32 = P.m[2][1] ; r33 = P.m[2][2] ;
% $$$ 
% $$$    /*                            [ r11 r12 r13 ]               */
% $$$    /* at this point, the matrix  [ r21 r22 r23 ] is orthogonal */
% $$$    /*                            [ r31 r32 r33 ]               */
% $$$ 
% $$$    /* compute the determinant to determine if it is proper */
% $$$ 
% $$$    zd = r11*r22*r33-r11*r32*r23-r21*r12*r33
% $$$        +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;  /* should be -1 or 1 */
% $$$ 
% $$$    if( zd > 0 ){             /* proper */
% $$$      ASSIF(qfac,1.0) ;
% $$$    } else {                  /* improper ==> flip 3rd column */
% $$$      ASSIF(qfac,-1.0) ;
% $$$      r13 = -r13 ; r23 = -r23 ; r33 = -r33 ;
% $$$    }
% $$$ 
% $$$    /* now, compute quaternion parameters */
% $$$ 
% $$$    a = r11 + r22 + r33 + 1.0l ;
% $$$ 
% $$$    if( a > 0.5l ){                /* simplest case */
% $$$      a = 0.5l * sqrt(a) ;
% $$$      b = 0.25l * (r32-r23) / a ;
% $$$      c = 0.25l * (r13-r31) / a ;
% $$$      d = 0.25l * (r21-r12) / a ;
% $$$    } else {                       /* trickier case */
% $$$      xd = 1.0 + r11 - (r22+r33) ;  /* 4*b*b */
% $$$      yd = 1.0 + r22 - (r11+r33) ;  /* 4*c*c */
% $$$      zd = 1.0 + r33 - (r11+r22) ;  /* 4*d*d */
% $$$      if( xd > 1.0 ){
% $$$        b = 0.5l * sqrt(xd) ;
% $$$        c = 0.25l* (r12+r21) / b ;
% $$$        d = 0.25l* (r13+r31) / b ;
% $$$        a = 0.25l* (r32-r23) / b ;
% $$$      } else if( yd > 1.0 ){
% $$$        c = 0.5l * sqrt(yd) ;
% $$$        b = 0.25l* (r12+r21) / c ;
% $$$        d = 0.25l* (r23+r32) / c ;
% $$$        a = 0.25l* (r13-r31) / c ;
% $$$      } else {
% $$$        d = 0.5l * sqrt(zd) ;
% $$$        b = 0.25l* (r13+r31) / d ;
% $$$        c = 0.25l* (r23+r32) / d ;
% $$$        a = 0.25l* (r21-r12) / d ;
% $$$      }
% $$$      if( a < 0.0l ){ b=-b ; c=-c ; d=-d; a=-a; }
% $$$    }
% $$$ 
% $$$    ASSIF(qb,b) ; ASSIF(qc,c) ; ASSIF(qd,d) ;
% $$$    return ;
% $$$ }

function [Cov_Matrix V D]=CovarianceMatrix(Componente_E,Componente_N,Componente_Z)

% COVARIANCE MATRIX : Determine the covariance matrix of a three component
%                    seismic sensor. The function also returns a diagonal 
%                    matrix D of eigenvalues and a matrix V whose columns 
%                    are the corresponding right eigenvectors.
%
% Input arguments  : (Component_E, Component_N ,Component_Z). The number of
%                    samples in each  component  must  be  the same lenght.
%                    The ordering of the three components must be specified
%                    as it is shown (E,N,Z).
%
% Output arguments : (Cov_Matrix, V, D). 3x3 Covariance matrix, real  and 
%                    symmetric. D is a diagonal matrix which contains the
%                    corresponding eigenvalues from the covariance matrix.
%                    V is a matrix whose columns are the corresponding 
%                    right eigenvectors.



%Autovariance E component (S11)
n=length(Componente_E);
sum=0;
s11=0;
    for i=1:n
        sum=Componente_E(i)*Componente_E(i);
        s11=s11+sum;
        i=i+1;
    end

    s11=s11/n;

%Autovariance N component (S22)
sum=0;
s22=0;
    for i=1:n
        sum=Componente_N(i)*Componente_N(i);
        s22=s22+sum;
        i=i+1;
    end

    s22=s22/n;

%Autovariance Z component (S33)
sum=0;
s33=0;
    for i=1:n
        sum=Componente_Z(i)*Componente_Z(i);
        s33=s33+sum;
        i=i+1;
    end
    
    s33=s33/n;

% Cross variance  E and N components (S12)
sum=0;
s12=0;
    for i=1:n
        sum=Componente_E(i)*Componente_N(i);
        s12=s12+sum;
        i=i+1;
    end

    s12=s12/n;

% Cross variance E and Z components (S13)
sum=0;
s13=0;
    for i=1:n
        sum=Componente_E(i)*Componente_Z(i);
        s13=s13+sum;
        i=i+1;
    end

    s13=s13/n;

% Cross variance N and E components (S21)
sum=0;
s21=0;
    for i=1:n
        sum=Componente_N(i)*Componente_E(i);
        s21=s21+sum;
        i=i+1;
    end

    s21=s21/n;

% Cross variance N and Z components (S23)
sum=0;
s23=0;
    for i=1:n
        sum=Componente_N(i)*Componente_Z(i);
        s23=s23+sum;
        i=i+1;
    end

    s23=s23/n;


% Cross variance Z and E components (S31)
sum=0;
s31=0;
    for i=1:n
        sum=Componente_Z(i)*Componente_E(i);
        s31=s31+sum;
        i=i+1;
    end

    s31=s31/n;

% Cross variance Z and N components (S32)
sum=0;
s32=0;
    for i=1:n
        sum=Componente_Z(i)*Componente_N(i);
        s32=s32+sum;
        i=i+1;
    end

    s32=s32/n;


Cov_Matrix=[s11 s12 s13;
            s21 s22 s23;
            s31 s32 s33];
        
% [V,D] = eig(A) returns diagonal matrix D of eigenvalues and matrix V whose
% columns are the corresponding right eigenvectors, so that A*V = V*D.

[V,D]=eig(Cov_Matrix);


end

% University of Granada - Final project of the Telecommunication engineering 
% degree - Signal Theory, Telematics and Communications Department (TSTC).
% Student : Reda Marzok Ben Omar.




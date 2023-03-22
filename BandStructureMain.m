%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                         %
%  This is the main file for Band Structure Plots. The folder which has this main file should have the    %
%  following accompanying files in it. Otherwise program will not work.                                   %
%                                                                                                         %
%    1) mystartdefaults.m                                                                                 %
%    2) bzpath.dat                                                                                        %
%    3) plot3Ddata.m                                                                                      %
%                                                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


source('mystartdefaults.m');  % Applies default setting from mystartdefaults.m file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Step 1 : Inputs

setpot = 1;                   % Choose one number from 1-14 to select the material from list below
nband=16;                     % No of bands to be stored in output file
cutoff = 21;                  % deal with only |G|^2 < cutoff (2*pi/spacing units)^2 is Hamiltonian
Gs_max = 11;                  % |G|^2  of highest non zero fourier coefficients in expanding potential




%           List of Materials

ListOfMaterials={"Si","Ge","Sn","GaP","GaAs","AlSb","InP","GaSb","InAs","InSb","ZnS","ZnSe","ZnTe","CdTe"};
MaterialName=char(ListOfMaterials(1,setpot));
fprintf("Selected material is %s.\n\n",MaterialName)

% Potentials from paper in ff (fourier forms) : First columns adjusts the top of filled band to zero
%          V0 VS3  VS8  VS11 VA3 VA4 VA11
ff(1,:) =   [-0.770437 -0.21 0.04 0.08 0 0 0];            % Si
ff(2,:) =   [-0.694179 -0.23 0.01 0.06 0 0 0];            % Ge
ff(3,:) =   [-0.500885 -0.20 0 0.04 0 0 0];               % Sn
ff(4,:) =   [-0.676246 -0.22 0.03 0.07 0.12 0.07 0.02];   % GaP
ff(5,:) =   [-0.651775 -0.23 0.01 0.06 0.07 0.05 0.01];   % GaAs
ff(6,:) =   [-0.509435 -0.21 0.02 0.06 0.06 0.04 0.02];   % AlSb
ff(7,:) =   [-0.561726 -0.23 0.01 0.06 0.07 0.05 0.01];   % InP
ff(8,:) =   [-0.51032 -0.22 0 0.05 0.06 0.05 0.01];       % GaSb
ff(9,:) =   [-0.523957 -0.22 0 0.05 0.08 0.05 0.03];      % InAs
ff(10,:) =  [-0.445297 -0.20 0 0.04 0.06 0.05 0.01];      % InSb
ff(11,:) =  [-0.466304 -0.22 0.03 0.07 0.24 0.14 0.04];   % ZnS
ff(12,:) =  [-0.448188 -0.23 0.01 0.06 0.18 0.12 0.03];   % ZnSe
ff(13,:) =  [-0.391512 -0.22 0.00 0.05 0.13 0.10 0.01];   % ZnTe
ff(14,:) =  [-0.309389 -0.20 0 0.04 0.15 0.09 0.04];      % CdTe

% lattice spacing of all the materials in angstrom
latticespacing = [5.43, 5.66, 6.49, 5.44, 5.64, 6.13, 5.86, 6.12, 6.04, 6.48, 5.41, 5.65, 6.07, 6.41];
spacing = latticespacing(1,setpot); % lattice spacing in angstrom for the selected material


%           Position of atoms in Primitive Cell

tau = zeros(3,2);
tau(:,1) = [0.125 0.125 0.125]';    % position of atom 1 in primitive cell
tau(:,2) = [-0.125 -0.125 -0.125]'; % position of atom 2 in primitive cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Step 2: Defining Unit Cell and Cell Volume in Cartesian Coordiante System

fprintf('FCC Lattice unit vectors in Cartesian Coordinates.\n')
a = zeros(3,3);
a(:,1) = [0.5 0.5 0.0]'; % direct lattice unit vector 1
a(:,2) = [0.0 0.5 0.5]'; % direct lattice unit vector 2
a(:,3) = [0.5 0.0 0.5]'; % direct lattice unit vector 3
printf("\n      a_1      a_2      a_3\n")
disp(a)

cell_volume = a(:,1)' * cross(a(:,2),a(:,3))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Step 3: Defining Reciprocal Lattice Vectors

g=zeros(4,3);

g(1:3,1) = cross(a(:,2),a(:,3))/cell_volume;
g(1:3,2) = cross(a(:,3),a(:,1))/cell_volume;
g(1:3,3) = cross(a(:,1),a(:,2))/cell_volume;

fprintf("\n\nFCC Reciprocal lattice unit vectors in Cartesian Coordinates.\n\n")
for i =1:3
  g(4,i) = g(1:3,i)' * g(1:3,i);
end
printf(" g_1 g_2 g_3\n")
disp(g)
min_norm = sqrt(min(g(4,:))) % minimal norm
nstep = floor(sqrt(cutoff)/min_norm) + 1;   %  Number of positive steps along each reciprocal lattice unit vector
printf("\n\nCutoff requires %d positive steps along each reciprocal lattice uniit vector.\n",nstep);

nodes = (2*nstep + 1)^3;  % Number of the reciprocal lattice vectors generated
printf("\n\nGenerate (2*%1d +1)^3 = %4d reciprocal lattice vectors\n\n",nstep,nodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Step 4: Generating reciprocal lattice vectors in all directions for calculations
G = zeros(5,nodes);
n=0;
for j = -nstep:nstep
  for k = -nstep:nstep
    for l = -nstep:nstep    % small L is iterating
      n++;
      G(1:3,n) = j*g(1:3,1)+k*g(1:3,2)+l*g(1:3,3);
      G(5,n) = G(1:3,n)'*G(1:3,n);
      G(4,n) = sqrt(G(5,n));
    end
  end
end

GT = sortrows(G',4); % sorting reciprocal lattice vectors by growing norm
G = GT';

kept = 1;
for n = 2:nodes
  if(G(5,n) <= cutoff)
    kept++;
  end
end

printf("%4d G vectors featuring |G|^2<cutoff\n",kept) % you should get 113 here
printf("  n             G(1)             G(2)             G(3)             |G|\n");
for i = 1:kept
printf("%3.6G %15.6G  %15.6G  %15.6G  %15.6G\n",i,G(1,i),G(2,i),G(3,i),G(4,i));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Step 5: Loading Brillouin Zone (bzpath.dat) Data File

datafile= 'bzpath.dat'; % 3D Brillouin zone path Data file name
specify_format = 'yes';
delim_in = ' ';
head_in = 1;

if (strcmp(specify_format,'yes'))
  [x, delim_out, head_out] = importdata(datafile, delim_in, head_in);
else
  [x,delim_out,head_out] = importdata(datafile);
end

sz = size(x.data);
nq = sz(1);       % No of q vectors
for iq = 1:nq
  is(iq) = x.data(iq,2);
  q(iq,1:5) = x.data(iq,3:7);
end

if head_out > 0
  x = x.data;
else
  x = x;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Step 6: For the kept reciprocal lattice vectors, calculate the value of V_G defined in Eq 4 (report)
%                                                                                                  Adopted from Monkhorst


if (setpot!=0) % fourier components of empty lattice potential
  ekinunit = 1.;
  printf("\nEnergy in (hbar^2/2*elm)*(2*pi/a)^2 units\n");
  for n=1:kept
    cvg(n) = 0+0i;
  end
end

if (setpot!=0) % fourier components of selected material pseudopotential
  printf("\nenergy in eV\n");
  ekinunit = ekin_fact*(2.d0*pi/spacing)^2 % kinectic energy in EV
  for n=1:kept
    sym = 0.;
    asym=0.;
    if(G(5,n) == 0)
      sym =ff(setpot,1)*Rydberg ; % here you can adjust zero of the potential to top of the valence band
      asym =0.;
    end
    if (G(5,n) == 3) % if (abs(G(5,n)-3)<tol)
      sym = ff(setpot,2)*Rydberg;
      asym = ff(setpot,5)*Rydberg;
    end
    if (G(5,n) == 4)
      sym =0.;
      asym = ff(setpot,6)*Rydberg;
    end
    if (G(5,n) == 8)
      sym = ff(setpot,3)*Rydberg;
      asym = 0.;
    end
    if (G(5,n) == 11)
      sym = ff(setpot,4)*Rydberg;
      asym =ff(setpot,7)*Rydberg;
    end
  argu = 2*pi * (G(1:3,n)'*tau(1:3,1));
  cvg(n) = cos(argu)*sym - 1i*sin(argu)*asym; % caution: sign of Im part
  end
end

printf("%4d G vectors featuring |G|^2<cutoff\n",kept) % you should get 113 here
  printf("  n             G(1)             G(2)             G(3)             |G|          |G|^2      Real Values   Imaginary Values\n");
for i = 1:kept
  printf("%3.6G %15.6G  %15.6G  %15.6G  %15.6G %15.6G %15.6G %15.6G \n",i,G(1,i),G(2,i),G(3,i),G(4,i),abs(G(4,i))^2,real(cvg(i)),imag(cvg(i)));
end

MyFileName=sprintf(['E',num2str(setpot),MaterialName,"-3D-EK-Diagram.dat"]); % Writing to a data files named with material
fprintf("\n\nName of the data file stored in computer for %s is: %s\n\n",MaterialName,MyFileName)
f2 = fopen(MyFileName,'w');

fprintf('\n Diagonalization loop over %4d wavevectors:',nq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Step 7: Initialize Hamiltonian matrix and assign potential energy values

G_diff = zeros(5,1);
H=zeros(kept,kept); % initialization of the Hamiltonian matrix
  %potential energy % you can place this potential energy part before the loop as optimization
  for j =1:kept
    for i =1:kept
      G_diff(1:3) = G(1:3,i) - G(1:3,j);
      G_diff(5) = G_diff(1:3)'*G_diff(1:3);
      if(G_diff(5) <= Gs_max)
      for k = 1:kept
        % if (abs(G_diff(1:3) - G(1:3,k))<[tol tol tol]')
        % would be better when comparing non-integer values
        if ((G_diff(1:3) - G(1:3,k)) == [0 0 0]')
          H(i,j) = cvg(k);
        end
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Step 8: Calculate the difference |k − G|^2 along the BZ exploration path and kinetic energy part of
%                   the hamiltonian matrix then diagonalize to find the eigenenergies

for iq =1:nq
  fprintf('%4d',iq); % This is counter to show on Command windows
  % Kinetic Energy
  for i = 1:kept
    for m=1:3
      p(m) = q(iq,m) - G(m,i);
    end
    H(i,i) = ekinunit * (p*p') + ff(setpot,1)*Rydberg; % + Adjusting top of the filled band to zero
  end

  % Hermiticity Check
  tol = 1e-10;
  if(!ishermitian(H,tol))
    printf("\nHamiltonian matrix not Hermitian : fatal error.\n");
    return;
  else
    % Diagonalization of Hamiltonian
    [v,ev]=eig(H) ;         % Diagonalization [eigenvectors, eigenvalues]
    E = real(diag(ev));     % Hermitian matrix features real eigenvalues
    [E,perm] = sort(E);     % Sorting eigenvalues in increasing order
    v = v(:,perm);          % Re-order eigenvectors in same order as permuted eigenvalues

    % Writing to file for band plotting
    bzs = 100 * (-1)^is(iq);
    if (kept < nband)
      nband = kept;
    end
    fprintf(f2, "%15.6G  %15.6G", q(iq,5), bzs);
    for i = 1:nband
      fprintf(f2,"%15.6G",E(i));
    end
    fprintf(f2, "\n");
    if (q(iq,4) == 0)
      for i=1:nband
        gamma(i) = E(i);
      end
    end
  end
end

% Writing the eigenenergies to a file
fclose(f2);
fprintf('\n\n%2d lowest eigenvalues at Gamma point: \n',nband);
for i = 1:nband
    fprintf('%4d %15.6G\n', i, gamma(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Step 9: Plot the eigenenergies against the BZ exploration path to obtain band structure
source('BS_Plot.m')

datafile= ['E',num2str(Setpot),'-',MaterialName,'-','MQ_r',num2str(MQ_r),'-','Smear',num2str(Smear),'-','DOSdata.dat'];


% Reading data file for specific material

specify_format= 'yes'; % If yes, the two next definitions will be used
delim_in= ' ';         % Expected column separator character
head_in =  1;          % Expected number of lines of header

labelx = 'Wavevector [2\pi/a units]'; % NB: Greek characters in TeX encoding
labely = 'Energy [eV]';

if(strcmp(specify_format,'yes'))
  [z,delim_out,head_out]=importdata(datafile,delim_in,head_in);
else
  [z,delim_out,head_out]=importdata(datafile);
end

fprintf('\nFOUND IN DATA FILE:\n');
fprintf('Column separator character = ''%s''\n',delim_out);
fprintf('Number of lines of header  = %d\n',head_out);
if(head_out>0)
  x=z.data; % If header found, file content read as "structure array"
            % (Display z to see what is meant)
else
  x=z;      % If header not found, file content read as numerical array
end

sz=size(x);
fprintf('abscissae = %d\n',sz(1))
fprintf('columns   = %d\n',sz(2))
x(sz(1,1),1)=0;
x(1,1)=0;
for i=1:length(x(:,1))
  if x(i,2)<0
    x1(i,1)=x(i,1);
    x1(i,2)=x(i,2);
  else
    x2(i,1)=x(i,1);
    x2(i,2)=x(i,2);
  endif
endfor
MaximumDensity=[max(x1(:,1)),max(x2(:,1))];
MaxDens=max(MaximumDensity);

myplot = figure()
##plot(x1(:,2),x1(:,1)); % Remember: Sourced plot defaults will apply!
h = area(x1(:,1)/MaxDens, x1(:,2));
set(h, 'FaceColor',[0.5 0.5 0.5]);
hold on
h1 = area(x2(:,1)/MaxDens, x2(:,2));
set(h1, 'FaceColor', 'w');
hold on
line(xlim,[0,0],"linestyle", "-.",'color','k','linewidth',0.2)
legend('Filled Bands','Empty Bands', "location", "southeast")
hold on


hold on;
ylim([-14,6]);

set(gca,'xtick',[]);
set(gca, 'ytick', [-14:2:6])
set(legend,'FontSize',6)
## mytitle=["Band Structure of ",MaterialName];
## title(mytitle,"interpreter","latex")
xlabel("Density of States (DOS) $[eV^{-1}]$","interpreter","latex")
ylabel("Energy (eV)","interpreter","latex")


			   %           Next lines must come after calling plot()


 				% OUTPUT FILE NAME

## dot=rindex(datafile, '.');      % Position of filetype extension delimiter
## root=substr(datafile,1,dot-1);  % Filename without filetype extension
## outfile=sprintf('%s.pdf',root);
##
## 				 %OUTPUT FILES
## print(outfile,'-dpdf');                % Basic pdf output
## print(outfile,'-dpdflatexstandalone'); % Combined pdf & LaTeX files
##
## fprintf('\nOUTPUT FILES: \n');
## fprintf('%s (basic pdf file)\n',outfile);
## fprintf('%s.tex and %s-inc.pdf (for perfect LaTeX processing)\n',root,root);

 waitfor(myplot); % Wait for closing graphical window

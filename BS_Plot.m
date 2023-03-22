datafile= ['E',num2str(setpot),MaterialName,"-3D-EK-Diagram.dat"]; % Reading data file for specific material

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
myplot = figure()
plot(x(:,1),x(:,2:end),'linewidth',1.5,'k'); % Remember: Sourced plot defaults will apply!
yline(0,"linestyle", "-.",'color','r','linewidth',0.5)
hold on;
ylim([-16,24]);

set(gca,'xtick',[0 0.5 1.5 1.75 2.5]);
set(gca, 'xticklabel',({"L","$\Gamma$","X","K","$\Gamma$"}));
set(gca, 'ytick', [-16:4:24])
##mytitle=["Band Structure of ",MaterialName];
##title(mytitle,"interpreter","latex")
xlabel("$k(2\pi/a)$","interpreter","latex")
ylabel("Energy (eV)","interpreter","latex")

			   %           Next lines must come after calling plot()


				% OUTPUT FILE NAME

dot=rindex(datafile, '.');      % Position of filetype extension delimiter
root=substr(datafile,1,dot-1);  % Filename without filetype extension
outfile=sprintf('%s.pdf',root);

				 %OUTPUT FILES
print(outfile,'-dpdf');                % Basic pdf output
print(outfile,'-dpdflatexstandalone'); % Combined pdf & LaTeX files

fprintf('\nOUTPUT FILES: \n');
fprintf('%s (basic pdf file)\n',outfile);
fprintf('%s.tex and %s-inc.pdf (for perfect LaTeX processing)\n',root,root);

waitfor(g); % Wait for closing graphical window

function DOUT=rainfallph(mat,tlim,OPT,nograph,dirspec)
% Rainwater Chemistry	
%
%       Plots Rainwater pH based on Ash Bucket codes.
%       Graphs are displayed under Environment -> Rainwater Chemistry on Webobs.
%       
%       Timesacles: '24hr', '07d', '30d', '1yr', 'all'
%       
%       This function was adapted from temp_loggers.m.
%
%       CGPS(MAT,TLIM,OPT,NOGRAPH) makes the following:
%           MAT = 1 (default) uses the local Matlab "past" data backup;
%           MAT = 0 forces to rebuilt the backup loading all the data from disk.
%           TLIM = DT or [T1;T2] plots a specific graph (extension '_xxx') for the 
%               last DT days, or on period betweendates T1 and T2, with vectorial format 
%               [YYYY MM DD] or [YYYY MM DD hh mm ss].
%           TLIM = 'all' plots a graphe of all the avalaible data time interval ('_all').
%           OPT.fmt = date format (see DATETICK).
%           OPT.mks = marker size.
%           OPT.cum = time cumulation interval for histograms (in days).
%           OPT.dec = decimation of data (in samples).
%           NOGRAPH = 1 (optionnal) do not plots any graphs.
%
%       DOUT = CGPS(...) returns a structure DOUT containing all the data from 
%       stations :
%           DOUT(i).code = station code (for station i)
%           DOUT(i).time = time vector (for station i)
%           DOUT(i).data = matrix of processed data (NaN = invalid data)
%
%
%   Authors: Arvid Ramdeane, MVO
%   Created: 2022-04-11
%   Updated: 

%diary rainfallph_log
%logfile = sprintf('');
%if nargin == 1
%	logfile = sprintf('rainfallph(%d)', mat);	
%end
%if nargin == 2
%	logfile = sprintf('rainfallph(%d,%s)', mat, tlim);
%end

%logfile_dir = '/var/log/webobs/matl

X = readconf;

if nargin < 1, mat = 1; end
if nargin < 2, tlim = []; end
if nargin < 3, OPT = 0; end
if nargin < 4, nograph = 0; end
if nargin < 5, dirspec = X.MKSPEC_PATH_WEB; end
if nargout > 0, nograph = 1; end

disp(nargout)

rcode = 'RAINFALLPH';
timelog(rcode,1);

% Initializing variables
G = readgr(rcode);
tnow = datevec(G.now);
samp = 10/1440;	% sampling rate (in days)
G.dsp = dirspec;
ST = readst(G.cod,G.obs);
ist = [find(~strcmp(ST.ali,'-') & ~strcmp(ST.dat,'-') & ~strcmp(ST.dat,''))];
aliases = [];

%datatype = '';	% name of the data directory (for now in 'pftp/<year>'
%folders - change if need be)
%samp = 1;	% sampling rate (in days)

pftp = sprintf('%s/%s',X.RACINE_FTP,G.ftp);
pdon = G.don;

%ydeb = 1998; %2003;
% f_tmp = '/tmp/.matlab-cgps-awk';
% f_tmp_header = '/tmp/.matlab-cgps-awk-header';

% f_tmp = sprintf('%s/.matlab-cgps-awk',pftp);
% f_tmp_header = sprintf('%s/.matlab-cgps-awk-header',pftp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding routine to read ash bucket data from one csv file
% Need only read the csv file once
datafile = sprintf('%s.csv', 'RainfallpH');
myfile = [pdon,'/',datafile];

% Check if file exists
% If file does not exist, exit program
if ~(exist(myfile, 'file'))
	disp('File does not exist. Exiting...');
	exit;
end

dd = importdata(myfile,',', 1); 

% Get number of columns in the file
% The first column will be Date, and the rest will be ash bucket codes
num_col = length(dd.textdata(1,:));
% Get dates in the date column and convert to datenum
% Need only do this once
t = datenum(dd.textdata(2:end,1), 'dd/mm/yyyy');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        



% Create a new directory to save generated images files

% Create new save directory
%{
new_dir = 'ashbucket_rainfallph';
current_save_dir = sprintf('%s/%s/%s',X.RACINE_FTP,G.ftp,X.MKGRAPH_PATH_FTP);
new_save_dir = sprintf('%s/%s',current_save_dir,new_dir);
% create new save directory if is doesnt not exists
if ~exist(new_save_dir, 'file')
	unix(sprintf('mkdir -p %s', new_save_dir));
end       
%}


for n = 1:length(ist)

	st = ist(n);
	scode = ST.cod{st};
	alias = ST.ali{st};
	sname = ST.nom{st};
    	data_file = ST.dat{st};%% column name in the csv file
	stitre = sprintf('%s - %s',alias, sname);
	stype = 'A';
    
	%disp(scode)

	%t = [];
	%d = [];
   
    	%% load the data
    
    	% Test: load Matlab backup file (if exists)
	%f_save = sprintf('%s/past/%s_past.mat',X.RACINE_OUTPUT_MATLAB,ts{n});
    	f_save = sprintf('%s/past/temp_%s_past.mat',X.RACINE_OUTPUT_MATLAB,alias);
	%f_save = 'C:\Users\Arvid\Documents\Adam\Webobs\rainwater_chemistry\temp_GALW_past.mat';
    	if mat && exist(f_save,'file')
		load(f_save,'t','d');
		fprintf('File: %s imported.\n',f_save)
	else
		disp('No Matlab backup. Loading all the data...');
        end
       	d = [];
       	for i=2:num_col
        	if strcmp(data_file, strtrim(dd.textdata{1,i}))
                	d = dd.data(:, i-1);
                	break
            	end
        end
        
		% if mat < 1 % load from text file

        	%datafile = 'RainfallpH_ABCS.csv';
        	%datafile = sprintf('%s.csv', data_file);
        
        	%myfile = [pftp,'\\',datafile];
        
       		%dd = csvread(file,'headerlines',1);
        	%dd = importdata(file);
        	%dd = importdata(myfile,','); 
        	%d = dd.data;
        	%t = datenum(dd.textdata(2:end,1));
        
    
        	%%% SAVE the data
        	save(sprintf('%s.part',f_save),'t','d');
	    	unix(sprintf('mv %s.part %s',f_save,f_save));
	    	disp(sprintf('File: %s updated.',f_save))
      
  	
 
   %%
    
	% Load the files (YYYY/SSSSYYYY.txt)
% 	for annee = ydeb:tnow(1)
% 		%f =
% 		%sprintf('%s/%s/%4d/%s%4d.txt',pftp,datatype,annee,ST.dat{st},annee
% 		%); 
%         
%         % removed datatype='final; folder. PJS, KP 2012-09-03
%         f = sprintf('%s/%4d/%s%4d.txt',pftp,annee,ST.dat{st},annee);
% 		% extracts only 12 column lines of the files
% 		if exist(f,'file') & ~unix(sprintf('awk ''{if (NF==11) print $0}'' %s > %s',f,f_tmp));
% 			dd = load(f_tmp);
% 			disp(sprintf('File: %s imported.',f));
% 			if ~isempty(dd)
% 				t = [t;datenum(dd(:,1:3))];
% 				d = [d;dd(:,5:10)];
% 			end
%         end
%         
%         % read header line to get component order for labels. PJS, KP
%         % 2012-08-29
%         if annee==tnow(1)
%         if exist(f,'file') & ~unix(sprintf('awk ''NR==1{print $6, $7, $8}'' %s > %s',f,f_tmp_header));
%         hdr{n}=textread(f_tmp_header,'%s');
%         end
%         end
%         
% 	end

% 	% Calibrates the data
% 	[d,C] = calib(t,d,ST.clb(st));
% 	nx = ST.clb(st).nx;
% 	so = 1:nx;
	nx = 1;

% 	% Exports in a single ASCII (for 'all' or 'xxx' only)
% 	if ischar(tlim)
% 		tt = datevec(t);
% 		f = sprintf('%s/%s.DAT',pftp,upper(scode));
% 		fid = fopen(f,'wt');
% 			fprintf(fid, '# DATE: %s\r\n', datestr(now));
% 			fprintf(fid, '# TITL: %s - %s\r\n',G.nom,stitre);
% 			fprintf(fid, '# SAMP: %d\r\n',round(samp*60*60*24));
% 			fprintf(fid, '# CHAN: YYYY MM DD');
% 			fmt = '%4d-%02d-%02d';
% 			for i = 1:nx
% 				fprintf(fid, ' %s_(%s)',C.nm{i},C.un{i});
% 				fmt = [fmt ' %0.3f'];
% 			end
% 			fprintf(fid,'\r\n');
% 			fmt = [fmt '\r\n'];
% 			fprintf(fid,fmt,[tt(:,1:6),d]');
% 		fclose(fid);
% 		disp(sprintf('File: %s updated.',f))
% 		clear tt
% 	end

	% Stores in other variables to prepare the synthesis graph
	eval(sprintf('d_%d=d;t_%d=t;',n,n));
	
	tlast(st) = t(end);

	% Interpr�tation des arguments d'entr�e de la fonction
	%	- t1 = temps min
	%	- t2 = temps max
	%	- structure G = param�tres de chaque graphe
	%		.ext = type de graphe (dur�e) "station_EXT.png"
	%		.lim = vecteur [tmin tmax]
	%		.fmt = num�ro format de date (fonction DATESTR) pour les XTick
	%		.cum = dur�e cumul�e pour les histogrammes (en jour)
	%		.mks = taille des points de donn�es (fonction PLOT)

    
	% Decoding argument TLIM
	if isempty(tlim)
        	ivg = 1:(length(G.ext)-2)
	end
	if ~isempty(tlim) & strcmp(tlim,'all')
		ivg = length(G.ext) - 1;
		G.lim{ivg}(1) = min(t);
    	end

	if ~isempty(tlim) & ~ischar(tlim)		
        if size(tlim,1) == 2		
            t1 = datenum(tlim(1,:));
			t2 = datenum(tlim(2,:));
        else		
            t2 = datenum(tnow);
			t1 = t2 - tlim;
        end
        
		ivg = length(G.ext);
		G.lim{ivg} = minmax([t1 t2]);
		if nargin > 2
			G.fmt{ivg} = OPT.fmt;
			G.mks{ivg} = OPT.mks;
			G.cum{ivg} = OPT.cum;
		end
		if ~nograph
			%f = sprintf('%s/%s/%s_xxx.txt',X.RACINE_WEB,dirspec,scode);
			%k = find(t>=G.lim{ivg}(1) & t<=G.lim{ivg}(2));
			%tt = datevec(t(k));
			%fid = fopen(f,'wt');
			%fprintf(fid, '# DATE: %s\r\n', datestr(now));
			%fprintf(fid, '# TITL: %s: %s\r\n',alias,sname);
			%fprintf(fid, '# SAMP: %d\r\n',samp*86400);
			%fprintf(fid, '# CHAN: YYYY MM DD HH NN %s_(%s) %s_(%s) %s_(%s)\r\n',C.nm{1},C.un{1},C.nm{2},C.un{2},C.nm{3},C.un{3});
			%fprintf(fid, '%4d-%02d-%02d %02d:%02d %0.2f %0.2f %0.2f\r\n',[tt(:,1:5),d(k,:)]');
			%fclose(fid);
			%disp(sprintf('Fichier: %s cr��.',f))
		end
	end
	% Send data to DOUT, for time period G.limit. -- Renvoi des donn�es dans DOUT, sur la p�riode de temps G.lim(end)
	if nargout > 0
		k = find(t>=G.lim{ivg(1)}(1) & t<=G.lim{ivg(1)}(2));
		DOUT(1).code = scode;
		DOUT(1).time = t(k);
		DOUT(1).data = d(k,:);
        DOUT(1).chan = C.nm;
		DOUT(1).unit = C.un;
	end

	% Don't draw graph is nograph==1 -- Si nograph==1, quitte la routine sans production de graphes
	if nograph == 1, ivg = []; end


	% ===================== Drawing of graphs ------------ Trac� des graphes
    
   % disp(ivg)
    
   for ig = ivg
      
        figure
		k = find(t>=G.lim{ig}(1) & t<=G.lim{ig}(2));
		if isempty(k)
			ke = [];
		else
			ke = k(end);
		end

            
			tk = t(k);
			dk = d(k);
        
% 		figure, clf, orient tall
% 		k = find(t>=G.lim{ig}(1) & t<=G.lim{ig}(2));
% 		if isempty(k)
% 			ke = [];
% 		else
% 			ke = k(end);
% 		end
% 
		% Etat de la station
                kacq = find(t>=G.lim{ig}(1) & t<=G.lim{ig}(2));
                if isempty(kacq)
			acqui = 0;
		else
			acqui = round(100*length(kacq)*samp/abs(t(kacq(end))-G.lim{ig}(1)));
		end
		%{
        if t(ke) >= G.lim{ig}(2)-G.lst
			etat = 0;
			for i = 1:nx
				if ~isnan(d(ke,i))
					etat = etat+1;
				end
			end
			etat = 100*etat/nx;
		else
			etat = 0;
		end
        %}
        etat = 0;
		% Title/info -------------------- Titre et informations
 		G.tit = gtitle(stitre,G.ext{ig});
        G.eta = [G.lim{ig}(2),etat,acqui];
 		etats(st) = etat;
 		acquis(st) = acqui;
		G.inf = {''};
        
 		if ig == 1 | strcmp(tlim,'all')
 			if ig == 1
				sd = '';
				for i = 1:nx
					%sd = [sd sprintf(', %1.1f %s', d(end,i),C.un{i})];
                    sd = [sd sprintf(', %1.1f', d(end,i))];
				end
				mketat(etat,tlast(st),sd(3:end),lower(scode),G.utc,acqui)
			end
		end

% 		if ~isempty(k)
% 			k1 = k(1);
% 			G.inf = {sprintf('Last measurement: {\\bf%s} {\\it%+d}',datestr(t(ke)),G.utc),' (min|moy|max)','','', ...
% 				sprintf('1. %s = {\\bf%+1.3f %s} (%+1.3f | %+1.3f | %+1.3f)',C.nm{1},d(ke,1)-d(k1,1),C.un{1},rmin(d(k,1)-d(k1,1)),rmean(d(k,1)-d(k1,1)),rmax(d(k,1)-d(k1,1))) ...
% 				sprintf('2. %s = {\\bf%+1.3f %s} (%+1.3f | %+1.3f | %+1.3f)',C.nm{2},d(ke,2)-d(k1,2),C.un{2},rmin(d(k,2)-d(k1,2)),rmean(d(k,2)-d(k1,2)),rmax(d(k,2)-d(k1,2))), ...
% 				sprintf('3. %s = {\\bf%+1.3f %s} (%+1.3f | %+1.3f | %+1.3f)',C.nm{3},d(ke,3)-d(k1,3),C.un{3},rmin(d(k,3)-d(k1,3)),rmean(d(k,3)-d(k1,3)),rmax(d(k,3)-d(k1,3))), ...
% 			};
% 		else
% 			k1 = [];
%         end
        
        if ~isempty(k)
			G.inf = {sprintf('Last measurement: {\\bf%s}',datestr(t(ke)))};%,' (min|mean|max)','','', ...
				%sprintf('Counts = {\\bf%+1.0f count} (%+1.0f | %+1.0f | %+1.0f)',d(ke),rmin(d(k)),rmean(d(k)),rmax(d(k))), ...
                %sprintf('Total %s events in last %s = %d',tp2d,T.nam{kk},sum(d(k)))
			%};
		else
			G.inf = '';
        end
		
		% loop for Relative Radial extension, Tangential displacement (clockwise) and Vertical displacements with error bars (in m)
		%for i = 1:3
			%subplot(6,1,(i-1)*2+(1:2)), extaxes %HERE NUMBER OF SUBPLOTS/STATIONS PLOTTED. 
% 			plot(repmat(t(k),[1,2])',(repmat(d(k,i)-d(k1,i),[1,2])+d(k,i+3)*[-1,1])','-','LineWidth',.1,'Color',.6*[1,1,1])
% 			hold on
% 			plot(t(k),d(k,i) - d(k1,i),'.','MarkerSize',G.mks{ig})
% 			hold off
            
            %plot(tk,dk,'.','MarkerSize',G.mks{ig})
            plot(tk,dk,'.','MarkerSize',10);
        
          
            set(gca,'XLim',[G.lim{ig}(1) G.lim{ig}(2)],'FontSize',8);
            set(gca, 'box', 'off');
 
			%datetick2('x',G.fmt{ig},'keeplimits')
            datetick2('x',1)
    


			%ylabel(sprintf('Relative %s (%s)',C.nm{i},C.un{i}))
            
            %THERE IS A PROBLEM, FOR NOW UNSOLVED, WITH the creation of hdr{5} and above (ie HERM and later). This is possibly due to the fact that HERM has no data 
            %at all for 2016. TO FIX THE PB TEMPORARILY, n IS FIXED TO 1 (corresponding to AIRS) SINCE IT ONLY IS USED HERE TO REFER TO EACH COMPONENT. 
            %SAME THING FOR ylabel command a few lines BELOW, KP, 08 oct 2016
            %switch hdr{n}{i}(1)
%             switch hdr{1}{i}(1)
%                 case 'R'
%                     cmp='Radial';
%                 case 'T'
%                     cmp='Tangential';
%                 case 'V'
%                     cmp='Vertical';
%             end
           ylabel('pH')
        %ylabel(sprintf('%s displacement %s',cmp,hdr{n}{i}(2:4))) %SEE PREVIOUS COMMENT
          %  ylabel(sprintf('%s displacement %s',cmp,hdr{1}{i}(2:4)))
            
			if isempty(find(~isnan(d(k)), 1)), nodata(G.lim{ig}), end
		%end

		tlabel(G.lim{ig},G.utc)
    
		ploterup
		filename_ext = 'rainfallph';
		image_filename = sprintf('%s_%s_%s',lower(scode),filename_ext,G.ext{ig});
		mkgraph(sprintf('%s_%s_%s',lower(scode),filename_ext,G.ext{ig}),G,OPT)
		
		% move generated image to new directory
		%unix(sprintf('mv %s/%s.* %s', current_save_dir, image_filename, new_save_dir)); 		

		close
	end
end

% 
% % ====================================================================================================
% % Graphs for all the network
% 
% stitre = sprintf('%s Network - final',G.nom);
% etat = mean(etats);
% acqui = mean(acquis);
% 
% for ig = ivg
% 
% 	figure, clf, orient tall
% 
% 	G.tit = gtitle(stitre,G.ext{ig});
% 	G.eta = [G.lim{ig}(2),etat,acqui];
% 	G.inf = {''};
% 	
% 	% time series (with light filter on the data)
% 	for i = 1:3
% 		subplot(6,1,(i-1)*2+(1:2)), extaxes
% 		hold on
% 		aliases = [];
% 		for n = 1:length(ist)
% 			eval(sprintf('d=d_%d;t=t_%d;',n,n));
% 			k = find(t>=G.lim{ig}(1) & t<=G.lim{ig}(2));
% 			if ~isempty(k)
% 				if length(k) > 1
% 					% filter the data (10*RMS on the derivate)
% 					drms = rstd(diff(d(k,i)));
% 					%disp(sprintf('%s: drms = %g',ST.ali{ist(n)},drms))
% 					kk = find(abs(diff(d(k,i)))>4*drms) + 1;
% 					d(kk,i) = NaN;
% 				end
% 				k1 = k(1);
% 				plot(t(k),d(k,i)-d(k1,i),'.','Color',scolor(n),'MarkerSize',G.mks{ig})
% 				aliases = [aliases,ST.ali(ist(n))];
% 			end
% 		end
% 		hold off
% 		set(gca,'XLim',[G.lim{ig}(1) G.lim{ig}(2)],'FontSize',8)
% 		box on
% 		datetick2('x',G.fmt{ig},'keeplimits')
%         
%         
%             switch hdr{n}{i}(1)
%                 case 'R'
%                     cmp='Radial';
%                 case 'T'
%                     cmp='Tangential';
%                 case 'V'
%                     cmp='Vertical';
%             end
%             ylabel(sprintf('%s displacement %s',cmp,hdr{n}{i}(2:4)))
%         
%         %ylabel(sprintf('Relative %s (%s)',C.nm{i},C.un{i}))
% 		
%         if length(find(~isnan(d(k,i))))==0, nodata(G.lim{ig}), end
% 		
% 		% legend
% 		xlim = get(gca,'XLim');
% 		ylim = get(gca,'YLim');
% 		nn = length(aliases);
% 		for n = 1:nn
% 			text(xlim(1)+n*diff(xlim)/(nn+1),ylim(2),aliases(n),'Color',scolor(n), ...
% 				'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',7,'FontWeight','bold')
% 		end
% 	end
% 
% 	tlabel(G.lim{ig},G.utc)
% 	ploterup
%     
% 	mkgraph(sprintf('%s_%s',rcode,G.ext{ig}),G,OPT)
% 	close
% 
% end
% 
% if isempty(tlim)
% 	mketat(etat,max(tlast),sprintf('%s %d stations',stype,length(etats)),rcode,G.utc,acqui)
% 	G.sta = [{rcode};lower(ST.cod(ist))];
% 	G.ali = [{rcode};ST.ali(ist)];
% 	G.ext = G.ext(1:end-1); % graph list: all except 'xxx'
% 	htmgraph(G);
% end
% 



%if isempty(tlim) && ~nograph
	mketat(etat,max(tlast),sprintf('%s %d stations',stype,length(etats)),rcode,G.utc,acqui)
	G.sta = [lower(ST.cod(ist))];
	% remaining G.sta
	for i=1:length(G.sta)
		G.sta{i} = sprintf('%s_rainfallph', G.sta{i});
		%disp(G.sta{i})
	end
	G.ali = [ST.ali(ist)];
	G.ext = G.ext(1:end-1); % graph list: all except 'xxx'
	htmgraph(G);
%end

timelog(rcode,2)

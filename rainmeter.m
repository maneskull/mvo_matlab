function DOUT=rainmeter(mat,tlim,OPT,nograph,dirspec)
%RAINMETER	Plots all graphs from rainmeter network.
%       RAINMETER without any option loads recent data from disk and updates all graphs on WEBOBS.
%
%       RAINMETER(MAT,TLIM,OPT,NOGRAPH) makes the following:
%           MAT = 1 (default) uses the local Matlab "past" data backup then update with new data files;
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
%       DOUT = RAINMETER(...) returns a structure DOUT containing all the data from 
%       stations :
%           DOUT(i).code = station code (for station i)
%           DOUT(i).time = time vector (for station i)
%           DOUT(i).data = matrix of processed data (NaN = invalid data)
%
%
%   Authors: F. Beauducel / IPGP + H. Odbert / MVO
%   Created: 2010-07-17
%   Updated: 2010-10-13
%
%   UPDATES
%
%   Authors: Arvid Ramdeane, MVO
%   Created: 2022-04-21
%
%   Reads and plots rainfall data from weather stations. Currently only MVO Helipad and Lees Yard data are plotted.
%   Other weather stations to come later.
%
%   RCODE: WXSTATIONS
% 
%   Data in csv format located in: /mvo/acqui/Raingauges/
%   This folder is binded to the network drive mvofls3 mounted in: /mnt/mvofls3/DVG_Data/WEBOBS_DATA/

X = readconf;

if nargin < 1, mat = 1; end
if nargin < 2, tlim = []; end
if nargin < 3, OPT = 0; end
if nargin < 4, nograph = 0; end
if nargin < 5, dirspec = X.MKSPEC_PATH_WEB; end
if nargout > 0, nograph = 1; end

rcode = 'WXSTATIONS';
timelog(rcode,1);

% Initializing variables
G = readgr(rcode);
tnow = datevec(G.now);

G.dsp = dirspec;
samp = 1/1440;	% sampling rate (in days)
pftp = sprintf('%s/%s',X.RACINE_FTP,G.ftp);
pdon = G.don;
ydeb = 2010;

ST = readst(G.cod,G.obs,1,[ydeb,1,1]);
ist = [find(~strcmp(ST.ali,'-') & ~strcmp(ST.dat,'-') & ~strcmp(ST.dat,'') & ST.ope>0)];
aliases = [];

for n = 1:length(ist)

	st = ist(n);
	scode = ST.cod{st};
	alias = ST.ali{st};
	sname = ST.nom{st};
	stitre = sprintf('%s: %s',alias,sname);
	stype = 'A';

	% Load data files: past loaded filenames are stored in a cell array "files".
	% New files are detected using setdiff function

	t = [];
	d = [];
	files = [];

	% Test: load Matlab backup file (if exists)
	f_save = sprintf('%s/past/%s_past.mat',X.RACINE_OUTPUT_MATLAB,lower(scode));
	if mat & exist(f_save,'file')
		load(f_save,'t','d','files');
		disp(sprintf('File: %s imported.',f_save))
	else
		disp('No Matlab backup. Loading all the data...');
	end

	%D = dir(sprintf('%s/%s_RainData2010+/%s_*.DAT',pdon,ST.dat{st}(1:3),ST.dat{st}));
	D = dir(sprintf('%s/%s_RainData/%s_*.csv',pdon, ST.dat{st}(1:3),ST.dat{st}));
	
	% Check if any files exist in directory D
	% If D is empty, possibly drive is not mounted, so no data files to read. Exit program...
	if isempty(D)
        	disp('No Data files in directory. Exiting...');
        	exit;
	end

	datafiles = cellstr(cat(1,D.name));
	newfiles = setdiff(datafiles,files);
	
	for i = 1:numel(newfiles)

		%f = sprintf('%s/%s_RainData2010+/%s',pdon,ST.dat{st}(1:3),newfiles{i});
		
		f = sprintf('%s/%s_RainData/%s',pdon,ST.dat{st}(1:3),newfiles{i});
		
		% imports file with format "11/04/2010,19:32,1.4" allowing:
		% 	- space instead of coma
		% 	- 2-digit year
		if exist(f,'file')
			
			% Date formated as: 24-Mar-21
			% Use importdata to import csv file
			%dd = importdata(f, ',', 1);
			%tt = datenum(dd.textdata(2:end, 1));
			%d = dd.data;
			%t = [t;tt];
			
			%[dd,mm,yy,hh,nn,data] = textread(f,'%n/%n/%n%n:%n%n%*[^\n]','headerlines', 1, 'whitespace',',\b\t');
			[dd,mm,yy,data] = textread(f,'%n/%n/%n%n%*[^\n]','headerlines', 1, 'whitespace',',\b\t');
			files = [files;newfiles(i)];
			disp(sprintf('File: %s imported.',f));
			if ~isempty(dd)
				k = find(yy < 100);
				yy(k) = yy(k) + 2000;
				%tt = datenum(yy,mm,dd,hh,nn,0) - 1/24;	% adjust time GMT+1 to GMT
				%tt = datenum(yy,mm,dd,hh,nn,0);  % adjust time GMT+1 to GMT
				tt = datenum(yy,mm,dd);  % adjust time GMT+1 to GMT
				d = [d;data];
				t = [t;tt];
			end
			
		end
	end

	% Sorts data, checks for overlaps and save the backup file
	if ~isempty(newfiles)
		[t,i] = sort(t);
		d = d(i);
		% note: after sorting, overlaps become repeated data
		k = find(diff(t) == 0);
		if ~isempty(k)
			t(k+1) = [];
			d(k+1) = [];
			disp(sprintf('Found %d overlaps. Deleted.',length(k))); %EDIT: change sprint to sprintf, Henry@MVO 15oct2010
		end

		save(sprintf('%s.part',f_save),'t','d','files');
		unix(sprintf('mv %s.part %s',f_save,f_save));
		disp(sprintf('File: %s updated.',f_save))
	end

	% data treatment: calibration and interpolation (to built a continuous 1-min vector)
	[d,C] = calib(t,d,ST.clb);
	nx = ST.clb.nx;

	ti = t(1):samp:t(end);	% 1-minute time vector
	di = xcum(t,d,ti)';	% 1-minute rain (no gap)

	if ~isempty(t)
		tlast(st) = t(end);
	else
		tlast(st) = NaN;
	end

	% Interpreting input arguments of the function
	%	- t1 = min time interval
	%	- t2 = max time interval
	%	- structure G = parameters of plots
	%		.ext = graph extension name (corresponding to time interval) "station_EXT.png"
	%		.lim = time vector limits [tmin tmax]
	%		.fmt = date format (function DATETICK) for XTick
	%		.cum = cumulative duration for histograms (in days)
	%		.mks = marker size (function PLOT)

	% Decoding argument TLIM
	if isempty(tlim)
		ivg = 1:(length(G.ext)-2);
	end
	if ~isempty(tlim) & strcmp(tlim,'all') | (isempty(tlim) & nograph == 1)
		ivg = length(G.ext)-1;
		G.lim{ivg}(1) = min(t);
	end
	if (~isempty(tlim) & ~ischar(tlim))
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
	end

	% Returns data in structure DOUT, within time interval G.lim(end)
	if nargout > 0
		k = find(t>=G.lim{ivg(1)}(1) & t<=G.lim{ivg(1)}(2));
		DOUT(1).code = scode;
		DOUT(1).time = t(k);
		DOUT(1).data = d(k,:);
		DOUT(1).chan = C.nm;
		DOUT(1).unit = C.un;
	end

	% If nograph==1, quits routine without making graphs
	if nograph == 1, ivg = []; end


	% ===================== Making graphs

	for ig = ivg

		figure
		k = find(t>=G.lim{ig}(1) & t<=G.lim{ig}(2));
		if isempty(k)
			ke = [];
			acqui = 0;
			tc = [];
			dc = [];
		else
			ke = k(end);
			acqui = round(100*length(k)*samp/abs(t(ke)-G.lim{ig}(1)));
			%tc = (G.lim{ig}(1)+.5*G.cum{ig}):G.cum{ig}:(G.lim{ig}(2)-.5*G.cum{ig});
			%tc = G.lim{ig}(1):samp:G.lim{ig}(2);	% 1-minute time vector
			%dk = xcum(t(k),d(k),tc);		% 1-minute rain
			ki = find(ti>=G.lim{ig}(1) & ti<=G.lim{ig}(2));
			tc = ti(ki);
			dc = filter(ones(round(G.cum{ig}/samp),1),1,di(ki));	% continuous rain
		end


		% Station status
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
        

		% Title and information
		G.tit = gtitle(stitre,G.ext{ig});
		G.eta = [G.lim{ig}(2),etat,acqui];
		etats(st) = etat;
		acquis(st) = acqui;
		if ig == 1 | strcmp(tlim,'all')
			if ig == 1
				sd = '';
				for i = 1:nx
					if ~isempty(d)
						sd = [sd sprintf(', %1.1f %s', d(end,i),C.un{i})];
					end
				end
				mketat(etat,tlast(st),sd(3:end),lower(scode),G.utc,acqui)
			end
		end

		if ~isempty(k)
			% computes some statistics
			d_sum = rsum(d(k));
			[d_max,i_max] = max(dc);
			t_max = tc(i_max(end));
			G.inf = {sprintf('Last measurement: {\\bf%s} {\\it%+d}',datestr(t(ke)),G.utc),'','','', ...
				sprintf('%s cumul. = {\\bf%1.0f %s}',C.nm{1},d_sum,C.un{1}) ...
				sprintf('%s max. per %s = {\\bf%1.0f %s} (last on %s)',C.nm{1},day2str(G.cum{ig}),d_max,C.un{1},datestr(t_max(end))) ...
			};
		else
			G.inf = '';
		end
		
		subplot(10,1,2:9), extaxes
		if ~isempty(k)
			%[ax,h1,h2] = plotyy(tc,dc,t(k),rcumsum(d(k)),'bar','plot');
			[ax,h1,h2] = plotyy(tc,dc,t(k),rcumsum(d(k)),'area','plot');
			colormap([0,1,1;0,1,1]), grid on
			ylim = get(ax(1),'YLim');
			set(ax(1),'XLim',G.lim{ig},'FontSize',8)
			ylim = get(ax(2),'YLim');
			set(ax(2),'XLim',G.lim{ig},'FontSize',8,'XTick',[])
			set(h2,'LineWidth',1.5)
			ylabel(sprintf('%s per %s (%s)',C.nm{1},day2str(G.cum{ig}),C.un{1}))
		
			ylabel(ax(1), 'Rainfall per day (mm)'); % ylabel
            		ylabel(ax(2), 'Total Rainfall (mm)');% ylabel
		else
			set(gca,'XLim',G.lim{ig},'FontSize',8)
			nodata(G.lim{ig})
		end
		
		%datetick2('x',G.fmt{ig},'keeplimits')
		datetick2('x',1,'keeplimits')
		tlabel(G.lim{ig},G.utc)

		if ~isempty(k)
			ploterup(ax(2))
		else
			ploterup(gca)
		end

		% Change plot figure width so y axis and label shows properly.
        	% Plot will be a little wider than usual.
		%gcf_position = get(gcf,'Position');
        	%set(gcf, 'Position', [gcf_position(1) gcf_position(2) 1200 gcf_position(3)]);  
    		
		mkgraph(sprintf('%s_%s',lower(scode),G.ext{ig}),G,OPT)
		close
	end
end


%if isempty(tlim) & ~nograph
	mketat(etat,max(tlast),sprintf('%s %d stations',stype,length(etats)),rcode,G.utc,acqui)
	G.sta = [lower(ST.cod(ist))];
	G.ali = [ST.ali(ist)];
	G.ext = G.ext(1:end-1); % graph list: all except 'xxx'
	htmgraph(G);
%end

timelog(rcode,2)


function [icarejcomps EEG_ICA_edited] = DefineRemoveArtificialICAs_EmoKids(cfg,EEG)
% This program aims to isolate artificial ICAs, e.g., VEOG, HEOG, and focal component, from the EEG data.
% However, "No automated method can accurately isolate artifacts without supervision (Chaumon et al., 2015),"
% which I totally agree; Therefore, if you don't trust the following automatic pipeline please use SASICA or 
% Adjust or other toolboxes to visually inspect all the ICAs, e.g., "PlotAndCheckICAsWithAdjust.m";
if ~isfield(cfg, 'focalcomp'),        cfg.focalcomp.enable = 1;                       end
if ~isfield(cfg, 'EOGcorr'),          cfg.EOGcorr.enable   = 1;                       end
if ~isfield(cfg, 'autocorr'),         cfg.autocorr.enable  = 1;                       end
if ~isfield(cfg, 'opts'),             cfg.opts.noplot      = 1;                       end
if ~isfield(cfg, 'plotarg'),          cfg.plotarg          = 0;                       end; plotarg = cfg.plotarg;       
if ~isfield(cfg, 'method'),           cfg.method           = 'both';                  end; method  = cfg.method;

%% run the ICAs detection with SASICA and Adjust with the program modified by WX;
% EOG correlation test;
cfg.EOGcorr.corthreshH    = 'auto'; %'auto 4' % M +/- 3*SD (99.x%) 2.5*SD (99%); M +/- 2*SD (95%);
cfg.EOGcorr.Veogchannames = 'VEOG';
cfg.EOGcorr.corthreshV    = 'auto'; %'auto 4'
cfg.EOGcorr.Heogchannames = 'HEOG';

% Focal component test; 
cfg.focalcomp.focalICAout = 'auto' %'auto 3';

% Muscle movement (auto  correlation); %The default is not to enable this function;
cfg.autocorr.dropautocorr = 'auto' %'auto 3';

% ADJUST
cfg.ADJUST.enable = 1;
[EEG1 output] = eeg_SASICA(EEG,cfg);
%% Remove the bad components;
% Results from SASICA
focalcomps  = find(EEG1.reject.SASICA.icarejfocalcomp); 
eogcomps    = find(EEG1.reject.SASICA.icarejchancorr);

if  cfg.autocorr.enable  == 0;
    mclcomps = [];
else
    mclcomps    = find(EEG1.reject.SASICA.icarejautocorr);
    gdsf_adjust = find(EEG1.reject.SASICA.icaADJUST.GDSF > EEG1.reject.SASICA.icaADJUST.soglia_GDSF);
    mclcomps    = intersect(mclcomps,gdsf_adjust);
end
% all bad components in SASICA
rejcomps_sasica = union(focalcomps, eogcomps);
if cfg.autocorr.enable  ~= 0;
    rejcomps_sasica = union(rejcomps_sasica,mclcomps);
end


% all bad components in ADJUST;

rejcomps_adjust = find(EEG1.reject.SASICA.icarejADJUST);
disp(['ADJUST reports these components are arteficial: ']);
disp(rejcomps_adjust);

disp(['SASICA reports these components are arteficial: ']);
disp(rejcomps_sasica);
%art = EEG1.reject.SASICA.icaADJUST.art;
%horiz = EEG1.reject.SASICA.icaADJUST.horiz;
%vert = EEG1.reject.SASICA.icaADJUST.vert;
%blink = EEG1.reject.SASICA.icaADJUST.blink;
%disc = EEG1.reject.SASICA.icaADJUST.disc;
%rejcomps_adjust = unique([art,horiz,vert,blink,disc]);

% Remove bad components identified by SASICA and ADJSUT;
switch method
	case 'sasica'
		icarejcomps = rejcomps_sasica;
	case 'adjust'
		icarejcomps = rejcomps_adjust;
	case 'common'
		icarejcomps = intersect(rejcomps_sasica,rejcomps_adjust);
	case 'both'
		icarejcomps = union(rejcomps_sasica,rejcomps_adjust);
end
% Both SASICA and ADJUST save their results in the gcompreject field, which is the default place for bad comps in EEGLAB;
EEG1.reject.gcompreject=EEG1.reject.gcompreject.*0;  % clear ADJUST's results;           
EEG1.reject.gcompreject(icarejcomps) = 1;  % Substitute it with the final rejected comps for the dataset;
if isempty(icarejcomps)
else
    disp(['Remove ' num2str(length(icarejcomps)) ' components: ' num2str(icarejcomps)]);
end
EEG_ICA_edited = pop_subcomp( EEG1, icarejcomps, 0);
if plotarg;
   % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','on'); 
   % eeglab redraw;
    pop_eegplot(EEG1,1,0,0,'','dispchans',10,'spacing',100);
    set (gcf,'name','Plot of EEG data before ICA rejection');
    EEG_ICA_edited.reject = EEG1.reject;
    pop_eegplot(EEG_ICA_edited,1,0,0,'','dispchans',10,'spacing',100,'children',gcf);
    set (gcf,'name','Plot of EEG data after ICA rejection');
    disp(['These components are removed: ' num2str(icarejcomps)]);
    SASICA(EEG);
end
clear EEG1 


disp('*********************** Finished ***********************');


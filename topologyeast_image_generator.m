function topologyeast_image_generator(tablename,outdir)
close all


TheMainTableForPredictionsfinel = import_top_file(tablename, "sheet1");
% Import data - cells 
for io=1:length(TheMainTableForPredictionsfinel(:,1))

    Prediction_chartl(TheMainTableForPredictionsfinel(io,:),outdir);
    
end

function [] =Prediction_chartl(Data_strain,outdir)
%get info from the main DB
ORF_name=Data_strain{1};%ORF_name_Column A
Gene_name=Data_strain{2};%Gene_name_Column B
Seq_lenght=Data_strain{4};%Seq_lenght_Column D 
signal_p=Data_strain{16};%signal_peptide_lengh TargetP(SP) column P
Probability_signal_p=Data_strain{17};%Probability of TargetP (SP) column Q
Treshold_for_SP=Data_strain{18};%Treshhol for singnalP 4.1 No SP when the probability is below the threshold column R
Validation_for_Nterminal=Data_strain{20};%Experimental validation of N’ “IN” orientation  column T
Validation_for_MTS=Data_strain{23};%Experimental validation of MTS column W
Probability_target_MTS=Data_strain{24};%Probability TargetP (SP) column X
lenght_target_MTS=Data_strain{25};%Lenght TargetP (SP) column W
Probability_MitoFates=Data_strain{26};%Probability MTS column Z
Lenght_MitoFates=Data_strain{27};%Lenght MTS column AA
Localization_of_N=Data_strain{28};%Localization of the N' terminus GFP Column AB
Localization_of_C=Data_strain{29};%Localization of the C' terminus GFP Column AC
Description=Data_strain{30};%Description column AD
Dubious=Data_strain{31};%Information about the dubiousness of the strain column AE
SPspoctopus_problem=Data_strain{33};%strains with problem with the SPoctopusan so that we should add acomment AG
%Retrotransposon=Data_strain{30};%Information about the dubiousness of the strain column AD
Matrix_predictions={};%Rows=type of prediction, Coloumn=Location of the characters, when 1=M 2=S 3=i 4=o.
%retrive info about Problem with spoctopus 
SP_prediction=Data_strain{13};
%adjast the scale for the protein size 
if Seq_lenght<= 1000
    Size_of_squre=140;
elseif (1300>=Seq_lenght)&&(Seq_lenght>1000)
    Size_of_squre=120;
elseif (1500>=Seq_lenght)&&(Seq_lenght>1300)
    Size_of_squre=90;
elseif Seq_lenght>=1500
    Size_of_squre=40;
end
%retrive topology 
for k=1:9    
    curent_prediction=Data_strain{6+k};
    Matrix_predictions{k,1}=strfind(curent_prediction,'M');
    Matrix_predictions{k,2}=strfind(curent_prediction,'S');
    Matrix_predictions{k,3}=strfind(curent_prediction,'i');
    Matrix_predictions{k,4}=strfind(curent_prediction,'o');
end
%build the figure
f=figure('Position',[310 0 1248 1150],'visible','off');
name_of_predictions={'TMHMM','TOPCONS','OCTOPUS','Philius','PolyPhobius','SCAMPI','SPOCTOPUS','HMMtop','memsat svm'};
cordination_for_y_axis={[0.904998429464832 0.819425986957506 0.0385416656732559 0.0319905207267305],...
                        [0.904998429464832 0.762518929106652 0.045833332122614 0.0319905207267305],...
                        [0.906560929464832 0.698532941528053 0.045833332122614 0.0319905207267305],...
                        [0.906560929464832 0.6280063462001 0.0343749991307656 0.0319905207267305],...
                        [0.906040096131499 0.568469516149344 0.0536458318897833 0.0319905207267305],...
                        [0.907602596131499 0.499934125276093 0.0385416656732559 0.0319905207267305],...
                        [0.907602596131499 0.440628747434641 0.0385416656732559 0.0319905207267305],...
                        [0.907602596131499 0.37270463657217 0.0385416656732559 0.0319905207267305],...
                       [0.905519262798165 0.307205555625162 0.0552083318432173 0.0319905207267305],};

%plot predection from the 9 Predictions programs                    
for l=1:9
    plot_seq=subplot(13,1,l+1);
    M=Matrix_predictions{l,1};
    S=Matrix_predictions{l,2};
    io=Matrix_predictions{l,3};
    o=Matrix_predictions{l,4};
    %soluble prediction
    scatter(M,ones(1,length(M))*1,Size_of_squre,[0.560784339904785 0.819607853889465 0.603921592235565],'s','filled')
    hold on
    scatter(S,ones(1,length(S))*1,Size_of_squre,[0.729411780834198 0.831372559070587 0.95686274766922],'s','filled')
    scatter(io,ones(1,length(io))*1,Size_of_squre,[0.929411768913269 0.686274528503418 0.686274528503418],'s','filled')
    scatter(o,ones(1,length(o))*1,Size_of_squre,[0.603921592235565 0.603921592235565 0.792156875133514],'s','filled')
    if length(o)==Seq_lenght || length(io)==Seq_lenght
        scatter(1:Seq_lenght,ones(1,Seq_lenght)*1,Size_of_squre,[0 0.800000011920929 0.800000011920929],'s','filled')
    end
    if Seq_lenght<=95%When the lenght of the sequence is too short under 70 BS we try to fix it by adding close squares.
        scatter(M-0.3,ones(1,length(M))*1,Size_of_squre,[0.560784339904785 0.819607853889465 0.603921592235565],'s','filled')
        scatter(S-0.3,ones(1,length(S))*1,Size_of_squre,[0.729411780834198 0.831372559070587 0.95686274766922],'s','filled')
        scatter(io-0.3,ones(1,length(io))*1,Size_of_squre,[0.929411768913269 0.686274528503418 0.686274528503418],'s','filled')
        scatter(o-0.3,ones(1,length(o))*1,Size_of_squre,[0.603921592235565 0.603921592235565 0.792156875133514],'s','filled')
       
        scatter(M+0.3,ones(1,length(M))*1,Size_of_squre,[0.560784339904785 0.819607853889465 0.603921592235565],'s','filled')
        scatter(S+0.3,ones(1,length(S))*1,Size_of_squre,[0.729411780834198 0.831372559070587 0.95686274766922],'s','filled')
        scatter(io+0.3,ones(1,length(io))*1,Size_of_squre,[0.929411768913269 0.686274528503418 0.686274528503418],'s','filled')
        scatter(o+0.3,ones(1,length(o))*1,Size_of_squre,[0.603921592235565 0.603921592235565 0.792156875133514],'s','filled')
        if length(o)==Seq_lenght || length(io)==Seq_lenght
        scatter((1:Seq_lenght)-0.3,ones(1,Seq_lenght)*1,Size_of_squre,[0 0.800000011920929 0.800000011920929],'s','filled')
        scatter((1:Seq_lenght)+0.3,ones(1,Seq_lenght)*1,Size_of_squre,[0 0.800000011920929 0.800000011920929],'s','filled')
        end
    end   
    Change_points=[];
    for in_sqe=2:Seq_lenght
        if Data_strain{6+l}(in_sqe) ~= Data_strain{6+l}(in_sqe-1)
            Change_points=[Change_points,in_sqe];
        end
    end
    %Problematic strains with SPoctopus
    if l==7 & Change_points > 1
        for CP=1:length(Change_points)-1
            if SP_prediction(Change_points(CP)-1)==SP_prediction(Change_points(CP+1)) & SP_prediction(Change_points(CP)-1)=='o'
                    annotation('textbox',...
                    [0.130608797290555 0.002739336492891 0.458662036042778 0.0725169673365439],...
                    'String',{'** Please note that the topology prediction generated by the TOPCONS batch runner of SPOCTOPUS for this gene is not logical and should be treated with caution. '},...
                    'FitBoxToText','off','EdgeColor','none');
                    %annotation('textbox',cordination_for_y_axis{7},'String','**SPOCTOPUS','FitBoxToText','on','EdgeColor','none','FontWeight','bold');
                    name_of_predictions{7}='**SPOCTOPUS';
            end
        end
    end
    Change_points=Change_points;% no idea way
    annotation('textbox',cordination_for_y_axis{l},'String',name_of_predictions{l},'FitBoxToText','on','EdgeColor','none','FontWeight','bold');
    box on
    range_seq=[-Seq_lenght/100 Seq_lenght+Seq_lenght/100];
    if Seq_lenght<=70
       xlim(range_seq+0.5);
    else
        xlim(range_seq);
    end
    ylim([0,1.5])
    if l==9
        set(gca,'xgrid','off');
    else
        set(gca,'xtick',[]);
    end
    set(gca,'ytick',[]);
    if not(isempty(Change_points))
        stem(Change_points-1,ones(1,length(Change_points))*1,'Color','r','Marker','o','LineStyle','-.',...
            'LineWidth',1,'MarkerSize',6,'MarkerFaceColor',[0.99 0.92 0.8],'MarkerEdgeColor','k')
        state=0.2;
        for in_sqe=1:length(Change_points)
            if in_sqe==1
                text_level=0.2;
                text(Change_points(in_sqe),text_level,num2str(Change_points(in_sqe)-1),'FontSize',8,'FontWeight','normal','Color','b');
                continue
            end
            if Change_points(in_sqe)-Change_points(in_sqe-1)<25 && state~=0.5
                text_level=0.5;
                text(Change_points(in_sqe),text_level,num2str(Change_points(in_sqe)-1),'FontSize',8,'FontWeight','normal','Color','b');
                state=0.5;
            else 
                text_level=0.2;
                text(Change_points(in_sqe),text_level,num2str(Change_points(in_sqe)-1),'FontSize',8,'FontWeight','normal','Color','b');
                state=0.2;
            end     
        end
    end
end

Original_tick=plot_seq.XTick;
if Seq_lenght<=Original_tick(end)        
    plot_seq.XTick(end)=Seq_lenght;
else
    plot_seq.XTick=[plot_seq.XTick,Seq_lenght];
end

if Seq_lenght >1300
    plot_seq.XTick(end-2)=[];
    plot_seq.XTick(end-1)=[];
elseif (plot_seq.XTick(end)-plot_seq.XTick(end-1))<30
    plot_seq.XTick(end-1)=[];
end

%plot signal_p data
   cur_plot=subplot(13,1,11);
   cur_pos=cur_plot.Position;
   set(gca,'Position',[cur_pos(1),cur_pos(2),cur_pos(3)-0.5,cur_pos(4)]);
if Probability_signal_p>=Treshold_for_SP
    if signal_p>=5
       scatter(1:signal_p,ones(1,signal_p)*1,Size_of_squre,[0.729411780834198 0.831372559070587 0.95686274766922],'s','filled')
       hold on
       stem(signal_p,1,'Color','r','Marker','o','LineStyle','-.',...
       'LineWidth',0.8,'MarkerSize',6,'MarkerFaceColor',[0.99 0.92 0.8],'MarkerEdgeColor','k')
       text(signal_p,0.2,num2str(signal_p),'FontWeight','normal','Color','b');
        %scatter(signal_p,ones(1,length(signal_p))*1,200,'y','s','filled')
    else
       Probability_signal_p='Below threshold';
    end
else 
    Probability_signal_p='Below threshold';
end
annotation('textbox',[0.4071875 0.210268240611297 0.223575211864407 0.0658057984536131],...
    'String',{['SignalP 4.1 (SP) probability: ' num2str(Probability_signal_p)]},...
    'EdgeColor','none','FontWeight','normal');
set(gca,'ytick',[])
xlim([-3 150]);
ylim([0,1.5]);
set(gca,'xtick',[]);
box on
%% probabilities target, MTS
subplot(13,1,12)
cur_plot=subplot(13,1,12);
cur_pos=cur_plot.Position;
set(gca,'Position',[cur_pos(1),cur_pos(2),cur_pos(3)-0.5,cur_pos(4)])
if lenght_target_MTS>=5
    scatter(1:lenght_target_MTS,ones(1,lenght_target_MTS)*1,Size_of_squre,[0.929411764705882 0.690196078431373 0.129411764705882],'s','filled')
    hold on
    stem(lenght_target_MTS,1,'Color','r','Marker','o','LineStyle','-.',...
    'LineWidth',0.8,'MarkerSize',6,'MarkerFaceColor',[0.99 0.92 0.8],'MarkerEdgeColor','k')
    text(lenght_target_MTS,0.2,num2str(lenght_target_MTS),'FontWeight','normal','Color','b');
else
    Probability_target_MTS='Below threshold';
end
annotation('textbox',...
    [0.4071875 0.143268240611297 0.223575211864407 0.0658057984536131],'String',{['TargetP (MTS) probability: ',num2str(Probability_target_MTS)]},...
    'EdgeColor','none','FontWeight','normal');
set(gca,'ytick',[])
xlim([-3 150]);
cur_plot.XTick=plot_seq.XTick;
ylim([0,1.5]);
set(gca,'xtick',[]);
box on
subplot(13,1,13)
cur_plot=subplot(13,1,13);
cur_pos=cur_plot.Position;
set(gca,'Position',[cur_pos(1),cur_pos(2),cur_pos(3)-0.5,cur_pos(4)])
if Probability_MitoFates>0.2 && Lenght_MitoFates > 5
    scatter(1:Lenght_MitoFates,ones(1,Lenght_MitoFates)*1,Size_of_squre,[0.929411764705882 0.690196078431373 0.129411764705882],'s','filled')
end
set(gca,'ytick',[])
xlim(range_seq);
cur_plot.XTick=plot_seq.XTick;
hold on
if Probability_MitoFates>0.2 && Lenght_MitoFates > 5
    stem(Lenght_MitoFates,1,'Color','r','Marker','o','LineStyle','-.',...
        'LineWidth',0.8,'MarkerSize',6,'MarkerFaceColor',[0.99 0.92 0.8],'MarkerEdgeColor','k')
    text(Lenght_MitoFates,0.2,num2str(Lenght_MitoFates),'FontWeight','normal','Color','b');
else Probability_MitoFates='Below threshold';
end
box on
xlim([-3 150]);
ylim([0,1.5]);
cur_plot.XTick=[0 50 100 150];
annotation('textbox',[0.4071875 0.083268240611297 0.223575211864407 0.0658057984536131],...
    'String',{['MitoFates (MTS) probability: ',num2str(Probability_MitoFates)]},...
    'EdgeColor','none','FontWeight','normal');
%plot legend of the Topology.
hold on
annotation('textbox',...
    [0.671041490480693 0.952699911615581 0.281486453408256 0.0504347351391072],'Color','k',...
    'String','Legend',...
    'FitBoxToText','off','EdgeColor','k',...
    'EdgeColor','none','FontWeight','bold','FontSize',12);
annotation('textbox',...
    [0.671041490480693 0.932699911615581 0.281486453408256 0.0504347351391072],'Color',[0.603921592235565 0.603921592235565 0.792156875133514],...
    'String','Out (Luminal or Extra-Cellular)',...
    'FitBoxToText','off','EdgeColor','k',...
    'EdgeColor','none','FontWeight','bold','FontSize',12);
annotation('textbox',...
    [0.671041490480693 0.912699911615581 0.281486453408256 0.0504347351391072],...
    'Color',[0.560784339904785 0.819607853889465 0.603921592235565],...
    'String','Trans Membrane Domain',...
    'FitBoxToText','off',...
    'EdgeColor','none','FontWeight','bold','FontSize',12);
annotation('textbox',[0.671041490480693 0.892699911615581 0.281486453408256 0.0504347351391072],...
    'Color',[0.729411780834198 0.831372559070587 0.95686274766922],...
    'String','Cleavable Signal peptide (SP)',...
    'FitBoxToText','off',...
    'EdgeColor','none','FontWeight','bold','FontSize',12);
annotation('textbox',...
    [0.671041490480693 0.872699911615581 0.151486453408256 0.0504347351391072],...
    'Color',[0.929411768913269 0.686274528503418 0.686274528503418],...
    'String','In (Cytosol)',...
    'FitBoxToText','off',...
    'EdgeColor','none','FontWeight','bold','FontSize',12);
annotation('textbox',...
    [0.671041490480693 0.852699911615581 0.281486453408256 0.0504347351391072],...
    'Color',[0.929411764705882 0.690196078431373 0.129411764705882],...
    'String','Mitochondrial Targeting Sequence (MTS) ',...
    'FitBoxToText','off',...
    'EdgeColor','none','FontWeight','bold','FontSize',12);
annotation('textbox',...
    [0.671041490480693 0.832699911615581 0.281486453408256 0.0504347351391072],...
    'Color',[0 0.800000011920929 0.800000011920929],...
    'String','Soluble',...
    'FitBoxToText','off',...
    'EdgeColor','none','FontWeight','bold','FontSize',12);
%print  the rest of data Name and description
%title
hold on
annotation('textbox',...
     [0.130608797290555 0.952699911615581 0.281370369376112 0.0496370656370712],...
    'String',{strcat(ORF_name,"  ",Gene_name)},...
    'EdgeColor','none','FontWeight','bold','FontSize',17);
annotation('textbox',...
    [0.130837745893887 0.953216374269006 0.13025 0.0198830409356725],...
    'String',{'SGD description:'},...
    'FitBoxToText','off',...
    'EdgeColor','none','FontWeight','bold');
annotation('textbox',...
    [0.130837745893887 0.877836257309942 0.531662254106113 0.0831627252762446],...
    'String',{Description},...
    'EdgeColor','none');
% expermintal validation Results:4-NO,3-vogtle 2009,2-Venna 2013,1-Yes both
%Experimental validation for n'
if Validation_for_MTS ==1
    annotation('textbox',[0.641041490480693 0.121806298553288 0.32646027505603 0.055868725868727],...
    'String',{'Experimental validation of MTS: \bfYes \rm, Venne et al 2013 , Vogtle et al 2009'},...
    'EdgeColor','none','FontSize',12);
elseif Validation_for_MTS==3
        annotation('textbox',[0.641041490480693 0.121806298553288 0.32646027505603 0.055868725868727],...
    'String',{'Experimental validation of MTS: \bfYes \rm, Vogtle et al 2009'},...
    'EdgeColor','none','FontSize',12);
elseif Validation_for_MTS==2
        annotation('textbox',[0.641041490480693 0.121806298553288 0.32646027505603 0.055868725868727],...
    'String',{'Experimental validation of MTS: \bfYes \rm, Venne et al 2013'},...
    'EdgeColor','none','FontSize',12);
elseif Validation_for_MTS==4
        annotation('textbox',[0.641041490480693 0.121806298553288 0.32646027505603 0.055868725868727],...
    'String',{'Experimental validation of MTS: \bfNo'},...
    'EdgeColor','none','FontSize',12);
end
    
if Validation_for_Nterminal==3
    annotation('textbox',...
    [0.641041490480693 0.0850943754029759 0.341777777777778 0.0551562876022501],'String',{"Experimental validation of N’ “in” orientation: \bfNot available"},...
    'EdgeColor','none','FontSize',12);
elseif Validation_for_Nterminal==2
    annotation('textbox',...
    [0.641041490480693 0.0850943754029759 0.341777777777778 0.0551562876022501],'String',{"Experimental validation of N’ “in” orientation: \bfYes"},...
    'EdgeColor','none','FontSize',12);
elseif Validation_for_Nterminal==1
        annotation('textbox',...
    [0.641041490480693 0.0850943754029759 0.341777777777778 0.0551562876022501],'String',{"Experimental validation of N’ “in” orientation: \bfNo"},...
    'EdgeColor','none','FontSize',12);
end
%print localiztion of GFp in N' and C'
annotation('textbox',...
    [0.641041490480693 0.228070175438597 0.0808495360427782 0.0500769577968333],...
    'String',{'Localization: '},...
    'EdgeColor','none','FontWeight','bold','FontSize',11);
%N' terminus
annotation('textbox',...
    [0.641041490480693 0.17906432748538 0.288229166666667 0.0526315789473685],...
    'String',{strcat("With Nt GFP tag: ",Localization_of_N)},...
    'EdgeColor','none','FontWeight','normal','FontSize',11);
%C' terminus
annotation('textbox',...
    [0.641041490480693 0.207000000000001 0.288229166666667 0.0478170744261607],...
    'String',{strcat("With Ct GFP tag: ",Localization_of_C)},...
    'EdgeColor','none','FontWeight','normal','FontSize',11);
ORF=Data_strain{1};
    Gene_name=Data_strain{2};
    Dubious=Data_strain{31};
    Retrotransopon=Data_strain{32};
    %save images to outdir
    if Dubious=='Y'
        saveas(gcf, strcat(outdir,filesep,ORF,"_",Gene_name,"(Dubious strain)",".png"));
    elseif Retrotransopon=='Y'        
        saveas(gcf, strcat(outdir,filesep,ORF,"_",Gene_name,"(Retrotransposon)",".png"));
    else
        saveas(gcf, strcat(outdir,filesep,ORF,"_",Gene_name,".png"));
    end
    close gcf
end
%import main DB
function TheMainTableForPredictionsfinel = import_top_file(workbookFile, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  THEMAINTABLEFORPREDICTIONSFINEL = IMPORTFILE(FILE) reads data from
%  the first worksheet in the Microsoft Excel spreadsheet file named
%  FILE.  Returns the data as a cell array.
%
%  THEMAINTABLEFORPREDICTIONSFINEL = IMPORTFILE(FILE, SHEET) reads from
%  the specified worksheet.
%
%  THEMAINTABLEFORPREDICTIONSFINEL = IMPORTFILE(FILE, SHEET, DATALINES)
%  reads from the specified worksheet for the specified row interval(s).
%  Specify DATALINES as a positive scalar integer or a N-by-2 array of
%  positive scalar integers for dis-contiguous row intervals.
%
%  Example:
%  TheMainTableForPredictionsfinel = importfile("D:\Predictions\The_Main_Table_For_Predictions_finel.xlsx", "Sheet1", [2, 6666]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 16-Feb-2022 15:11:36

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [2, 6666];
end

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 41);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "A" + dataLines(1, 1) + ":AO" + dataLines(1, 2);

% Specify column names and types
opts.VariableNames = ["ORF", "geneName", "MWInKd", "Length", "SequenceFromSGD", "FASTAFormat", "Tmhmm", "TOPCONS", "OCTOPUS", "Philius", "PolyPhobius", "SCAMPI", "SPOCTOPUS", "HMMtop", "memsat_svm", "signal_peptide_lenghSignal_P41SP", "Probability_of_signalPSignal_P41SP", "ThresholdForSignal_P41SP", "SPPredictionClassIsNotShown", "ValueForTheProgramOfNINOrientation3NotAvalible2Yes1No", "Venne2013IsNotShown", "Vogtle2009IsNotShown", "ResultsForExperimentalValidationForMTS4NOForBothOfThem3vogtle20", "ProbabilitytargetPMTS", "LengthTargetPMTS", "ProbabilityMitofatesMTS", "LengthMItoFatesMTS", "LocalizationByNGFP", "LocalizationByCGFP", "FullDescreptionBySGD", "DubiousYN", "RetrotransposonYN", "StrainsProblematicWithSpoctopus001ES1", "VarName34", "VarName35", "VarName36", "VarName37", "length1", "length2", "P", "VarName41"];
opts.VariableTypes = ["char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char"];

% Specify variable properties
opts = setvaropts(opts, ["ORF", "geneName", "MWInKd", "Length", "SequenceFromSGD", "FASTAFormat", "Tmhmm", "TOPCONS", "OCTOPUS", "Philius", "PolyPhobius", "SCAMPI", "SPOCTOPUS", "HMMtop", "memsat_svm", "signal_peptide_lenghSignal_P41SP", "Probability_of_signalPSignal_P41SP", "ThresholdForSignal_P41SP", "SPPredictionClassIsNotShown", "ValueForTheProgramOfNINOrientation3NotAvalible2Yes1No", "Venne2013IsNotShown", "Vogtle2009IsNotShown", "ResultsForExperimentalValidationForMTS4NOForBothOfThem3vogtle20", "ProbabilitytargetPMTS", "LengthTargetPMTS", "ProbabilityMitofatesMTS", "LengthMItoFatesMTS", "LocalizationByNGFP", "LocalizationByCGFP", "FullDescreptionBySGD", "DubiousYN", "RetrotransposonYN", "StrainsProblematicWithSpoctopus001ES1", "VarName34", "VarName35", "VarName36", "VarName37", "length1", "length2", "P", "VarName41"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["ORF", "geneName", "MWInKd", "Length", "SequenceFromSGD", "FASTAFormat", "Tmhmm", "TOPCONS", "OCTOPUS", "Philius", "PolyPhobius", "SCAMPI", "SPOCTOPUS", "HMMtop", "memsat_svm", "signal_peptide_lenghSignal_P41SP", "Probability_of_signalPSignal_P41SP", "ThresholdForSignal_P41SP", "SPPredictionClassIsNotShown", "ValueForTheProgramOfNINOrientation3NotAvalible2Yes1No", "Venne2013IsNotShown", "Vogtle2009IsNotShown", "ResultsForExperimentalValidationForMTS4NOForBothOfThem3vogtle20", "ProbabilitytargetPMTS", "LengthTargetPMTS", "ProbabilityMitofatesMTS", "LengthMItoFatesMTS", "LocalizationByNGFP", "LocalizationByCGFP", "FullDescreptionBySGD", "DubiousYN", "RetrotransposonYN", "StrainsProblematicWithSpoctopus001ES1", "VarName34", "VarName35", "VarName36", "VarName37", "length1", "length2", "P", "VarName41"], "EmptyFieldRule", "auto");

% Import the data
TheMainTableForPredictionsfinel = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "A" + dataLines(idx, 1) + ":AO" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    TheMainTableForPredictionsfinel = [TheMainTableForPredictionsfinel; tb]; %#ok<AGROW>
end

%% Convert to output type
TheMainTableForPredictionsfinel = table2cell(TheMainTableForPredictionsfinel);
numIdx = cellfun(@(x) ~isnan(str2double(x)), TheMainTableForPredictionsfinel);
TheMainTableForPredictionsfinel(numIdx) = cellfun(@(x) {str2double(x)}, TheMainTableForPredictionsfinel(numIdx));
end
end
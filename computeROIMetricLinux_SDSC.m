function computeROIMetricLinux_SDSC(reflist, tstlist, resultfile, winsuffix)
%computeROIMetricLinux('reflist.csv', 'autolist.csv', 'prostateTestMetrics3.csv')
% reflist = 'reflist.csv';
% tstlist = 'forseglist.csv';
% resultfile = 'test.csv';
if(nargin<4)
    winsuffix='';
end
reftable=readtable(reflist, 'delimiter', ',');
tsttable=readtable(tstlist, 'delimiter', ',');
fp=fopen(resultfile, 'w');
reftablesize=size(reftable);
tsttablesize=size(tsttable);
numnonstrs=4;
numstrs=tsttablesize(2)-numnonstrs;
dice=zeros(tsttablesize(1), numstrs);
fprintf(fp, 'tstfname, reffname ');

SurfaceDSCThreshold = 0.3;% 2mm threshold for surfaceDice

indlabel={'dice', 'recall', 'precision', 'hausdorff90','hausdorff95', 'hausdorff98', 'meanSurDist', ['SurfaceDice' num2str(SurfaceDSCThreshold)]};
inds=numel(indlabel); %dice, recall, precision, hausdorff
for i=1:numstrs
    stname=tsttable.Properties.VariableNames{numnonstrs+i};
    for j=1:inds
        fprintf(fp, ',%s', [stname, '_', indlabel{j}]);
    end
end
fprintf(fp, '\n');

for i=1:tsttablesize(1) 
    fname=tsttable{i,1}{1}
    
    found=0;
    
    for j=1:reftablesize(1)
        if(strcmp(reftable{j,2}{1}, tsttable{i,2}{1}))
            found=j;
            break;
        end
    end
    if(~found)
        continue;
    end
    
    reffname=reftable{found,1}{1};
    if(nargin>=4)
        infotst=dicominfo([winsuffix ], 'UseVRHeuristic', false);
        inforef=dicominfo([winsuffix reffname(6:end)], 'UseVRHeuristic', false);
    else
        infotst=dicominfo(fname, 'UseVRHeuristic', false);
        inforef=dicominfo(reffname, 'UseVRHeuristic', false);
    end
    %infotst=dicominfo(['Y:' fname(6:end)], 'UseVRHeuristic', false);
    %inforef=dicominfo(['Y:' reffname(6:end)]);
    fprintf(fp, '%s, %s ', fname, reffname);
    hd=zeros(numstrs, 11);
    for k=(numnonstrs+1):tsttablesize(2)
        stname=tsttable.Properties.VariableNames{k}
        structtst=[];
        structref=[];
        if(iscell(tsttable{i,k}) && iscell(reftable{found,k}))
        tstname=tsttable{i,k}{1};
        refname=reftable{found,k}{1};
        
        structtst=getStructbyName(infotst, tstname);
        structref=getStructbyName(inforef, refname);
        end
        if(isempty(structref)||isempty(fieldnames(structref)))
            inter(k-numnonstrs)=0;
            refonly(k-numnonstrs)=0;
            tstonly(k-numnonstrs)=1; 
            union(k-numnonstrs)=1;
            meanD(k-numnonstrs)=NaN;
            hd(k-numnonstrs,:)=NaN;
            surfaceDice(k-numnonstrs)=NaN;
            continue;
        end
        if(isempty(structtst)||isempty(fieldnames(structtst)))
            inter(k-numnonstrs)=0;
            refonly(k-numnonstrs)=1;
            tstonly(k-numnonstrs)=0; 
            union(k-numnonstrs)=1;
            meanD(k-numnonstrs)=NaN;
            hd(k-numnonstrs,:)=NaN;
            surfaceDice(k-numnonstrs)=NaN;
            continue;
        end
        if(strcmpi(stname, 'SpinalCord')||strcmpi(stname, 'Rectum')||strcmpi(stname, 'FemoralHead_R')||strcmpi(stname, 'FemoralHead_L')||strcmpi(stname, 'SpinalCanal')||strcmpi(stname, 'Aorta')||strcmpi(stname, 'V_Venacava_I'))
            flag='overlap'  % overlap
        elseif(strcmpi(stname, 'Esophagus'))
            flag='refonly'
        else
            flag='default'
        end
structref
structtst
        [inter(k-numnonstrs), refonly(k-numnonstrs), tstonly(k-numnonstrs), union(k-numnonstrs), hd(k-numnonstrs,:), meanD(k-numnonstrs), surfaceDice(k-numnonstrs)]=computeOverlap(structref, structtst, flag, SurfaceDSCThreshold);  % 2mm threshold for surfaceDice
        
    end
    dice(i,:)=2*inter./(inter*2+refonly+tstonly);
    tpr(i,:)=inter./(inter+refonly);
    ppv(i,:)=inter./(inter+tstonly);
    hd90s(i,:)=hd(:,1)';
    hd95s(i,:)=hd(:,6)';
    hd98s(i,:)=hd(:,9)';
    meands(i,:)=meanD;
    surfdices(i,:)=surfaceDice;
    for ii=1:numstrs
        fprintf(fp, ',%f,%f,%f,%f,%f,%f,%f,%f', dice(i,ii), tpr(i, ii), ppv(i, ii),hd90s(i,ii),hd95s(i,ii), hd98s(i,ii), meands(i,ii),surfdices(i,ii));
    end
    fprintf(fp, '\n');
end
fclose(fp);
end

function struct=getStructbyName(roi, name)
temp1=fieldnames(roi.StructureSetROISequence);
temp=fieldnames(roi.ROIContourSequence);
numstrsc=length(temp1);
numstrsc1=length(temp);
roiNames{numstrsc}=[];
roiNums(numstrsc)=0;
foundc=0;
struct=[];
for kc=1:numstrsc
     if(strcmpi(roi.StructureSetROISequence.(temp1{kc}).ROIName,name))
         foundc=1;
         foundroinum=roi.StructureSetROISequence.(temp1{kc}).ROINumber;
         break;
     end
end
if (foundc)
    for jc=1:numstrsc1
        if(isfield(roi.ROIContourSequence.(temp{jc}), 'ReferencedROINumber')&&(foundroinum==roi.ROIContourSequence.(temp{jc}).ReferencedROINumber)&&isfield(roi.ROIContourSequence.(temp{jc}), 'ContourSequence'))
            struct=roi.ROIContourSequence.(temp{jc}).ContourSequence;
            break
        end
    end
end
end


function [inter, refonly, tstonly, union, hd, meanD, surfaceDice]=computeOverlap(structref, structtst,flag, surfaceDiceThresholdinCM)


[refmins, refmaxs, refslicethickness, structref]=analyStruct(structref);
[tstmins, tstmaxs, tstslicethickness, structtst]=analyStruct(structtst);
if nargin<3
    flag='default';
end
allmins=min(refmins, tstmins); % find the allminimum value between refmins and tstmins in x,y,z, seperatively
allmaxs=max(refmaxs, tstmaxs);% find the allmaximum value between refmaxs and tstmaxs in x,y,z, seperatively
if(strcmpi(flag, 'refonly')) % if refonly
    allmins(3)=refmins(3);
    allmaxs(3)=refmaxs(3);
elseif(strcmpi(flag, 'overlap'))%if overlap 
        allmins(3)=max(refmins(3), tstmins(3));
        allmaxs(3)=min(refmaxs(3), tstmaxs(3));
end

header.x_pixdim=1/10; % x pixel dimension (1/10 means the unit of mm?)
header.y_pixdim=1/10; % y pixel dimension
header.z_pixdim=max([refslicethickness, tstslicethickness])/10; % z pixel dimension is the 1/10 * max between reference thickness and test thickness
% creater the header struct, with fields x_pixdim:0.1, y_pixdim:0.1,  z_pixdim: 0.2500

if(header.z_pixdim==0) % z pixel dimension cannot ==0, if it is, get it a value 0.1
    header.z_pixdim=0.1;
end
header.x_start=allmins(1)/10;
header.y_start=allmins(2)/10;
% x, y start from the allmins(1), allmins(2)
header.z_start=-allmaxs(3)/10;
% z start from the -allmaxs(3)/10?
header.x_dim=ceil((allmaxs(1)-allmins(1))/10/header.x_pixdim)+2;
header.y_dim=ceil((allmaxs(2)-allmins(2))/10/header.y_pixdim)+2;
header.z_dim=ceil((allmaxs(3)-allmins(3))/10/header.z_pixdim)+1;
% ceil():Round toward positive infinity
refmask=computeMaskfromStruct(structref, header, 'HFS'); % compute the mask from struct 
[Bref, ~]=bwboundaries2Dstackfast(refmask);
refpts=[];
tstpts=[];
for k=1:length(Bref)
    refpts=[refpts; Bref{k}];
end
tstmask=computeMaskfromStruct(structtst, header, 'HFS');
[Btst, ~]=bwboundaries2Dstackfast(tstmask);
if(isempty(Btst))
    hd=ones(1,11)*NaN;
    inter=NaN;
    union=NaN;
    tstonly=NaN;
    refonly=NaN;
    meanD=NaN;
    return;
end
for k=1:length(Btst)
    tstpts=[tstpts; Btst{k}];
end
refpts=refpts.*(ones(size(refpts,1),1)*[header.x_pixdim, header.y_pixdim, header.z_pixdim]);
tstpts=tstpts.*(ones(size(tstpts,1),1)*[header.x_pixdim, header.y_pixdim, header.z_pixdim]);
[hd, meanD,~, dp, dq]=HausdorffDistPctile(refpts, tstpts,[90:100]);

surfaceDice=(sum(dp<surfaceDiceThresholdinCM)+sum(dq<surfaceDiceThresholdinCM))/(numel(dp)+numel(dq));

hd95=hd(6);
hd98=hd(9);
bitMask=zeros([header.x_dim, header.y_dim, header.z_dim], 'uint8');
bitMask=bitset(bitMask, 1, refmask);
bitMask=bitset(bitMask, 2, tstmask);
inter=sum(bitMask(:)==3);
refonly=sum(bitMask(:)==1);
tstonly=sum(bitMask(:)==2);
union=sum(bitMask(:)>0);

end

function [mins, maxs, slicethickness, instruct]=analyStruct(instruct)
    temp=fieldnames(instruct);

    slicez(length(temp))=0;
    for i=1:length(temp)
        if(isempty(instruct.(temp{i})))
            instruct=rmfield(instruct, temp{i});
        else
        numpoints=instruct.(temp{i}).NumberOfContourPoints;
        pts=zeros(3,numpoints);
        pts(:)=instruct.(temp{i}).ContourData(:);
        if(numpoints>1)
            tmpmins=min(pts');
            tmpmaxs=max(pts');
        else
            tmpmins=pts';
            tmpmaxs=pts';
        end
        slicez(i)=tmpmins(3);
        if(i<2)
            mins=tmpmins;
            maxs=tmpmaxs;
        else
            mins=min(mins, tmpmins);
            maxs=max(maxs, tmpmaxs);
        end
        end
    end
    if(length(temp)>0&&isfield(instruct.Item_1, 'ContourSlabThickness'))
            slicethickness=instruct.(temp{i}).ContourSlabThickness;
    else if(~isfield(instruct.Item_1, 'ContourSlabThickness')&& length(temp)>1)
            slicethickness=median(diff(sort(unique(slicez), 'ascend')));
        else
            slicethickness=0;
        end
    end
end


function varargout = Pcttool(varargin)
% PCTTOOL M-file for Pcttool.fig
%      PCTTOOL, by itself, creates a new PCTTOOL or raises the existing
%      singleton*.
%
%      H = PCTTOOL returns the handle to a new PCTTOOL or the handle to
%      the existing singleton*.
%
%      PCTTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PCTTOOL.M with the given input arguments.
%
%      PCTTOOL('Property','Value',...) creates a new PCTTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Pcttool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Pcttool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Pcttool

% Last Modified by GUIDE v2.5 03-Jul-2012 11:22:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Pcttool_OpeningFcn, ...
                   'gui_OutputFcn',  @Pcttool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Pcttool is made visible.
function Pcttool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Pcttool (see VARARGIN)

% Choose default command line output for Pcttool
handles.output = hObject;


%My initialization


axes(handles.axes_image);
set(gca,'DataAspectRatio', [1 1 1]);
set(gcf,'WindowButtonMotionFcn', @Tactool_Cursorpos);

%Image fields
handles.activeData = 'none'; %Could also be 'raw','preprocessed',and 'residue'.
handles.activeImage = 'none'; %Could also be 'cbv','cbf','mtt','bbbp'
handles.isSequence = false; %True if current display is a sequence of images

handles.haveData = false;
handles.nImages = 0;
handles.currentImageNo = 1;
handles.imageData = [];
handles.mapData = [];
handles.image_height = 0;
handles.image_width = 0;

handles.havePreprocessed = false;
handles.haveResidueData = false;
handles.havePerfusion = false;

handles.cbv.haveData = false;
handles.cbv.data = [];
handles.cbf.haveData = false;
handles.cbf.data = [];
handles.mtt.haveData = false;
handles.mtt.data = [];
handles.ttp.haveData = false;
handles.ttp.data = [];
handles.bbbp.haveData = false;
handles.bbbp.data = [];
handles.patlak.haveData = false;
handles.patlak.xmap = [];
handles.patlak.ymap = [];
handles.patlak.r = [];
handles.patlak.bbbp = [];

%For the TAC/Patlak dump
handles.x = 0;
handles.y = 0;

%Initialize the settings struct with standard values
handles.settings.tissue_density = 1.05;
handles.settings.hematocrit = 0.73;
handles.settings.lower_threshold = 0;
handles.settings.upper_threshold = 120;
handles.settings.sfilter_size = 5;
handles.settings.sfilter_stddev = 0.5;
handles.settings.tfilter_size = 5;
handles.settings.tfilter_stddev = 0.5;
handles.settings.aif_x = 1; 
handles.settings.aif_y = 1;
handles.settings.vof_x = 1;
handles.settings.vof_y = 1;
handles.settings.truncation = 0.3;
handles.settings.block_m = 2;
handles.settings.first_frame = 2;
handles.settings.last_frame = -1;
handles.settings.tracer_factor = 1;
handles.settings.enable_sfilter = true;
handles.settings.enable_nldf = false;
handles.settings.enable_tfilter = true;
handles.settings.enable_segment = true;
handles.settings.enable_subtractbase = true;
handles.settings.enable_normalizeAIF = true;
handles.settings.enable_hematocrit = true;
handles.settings.min_preprocessed = 0;
handles.settings.max_preprocessed = Inf;
%Perfusion settings
handles.pctsettings.min_cbv = 0;
handles.pctsettings.max_cbv = Inf;
handles.pctsettings.min_cbf = 0;
handles.pctsettings.max_cbf = Inf;
handles.pctsettings.min_mtt = 0;
handles.pctsettings.max_mtt = 40.95;
handles.pctsettings.min_ttp = 0; %TODO: Change these values!
handles.pctsettings.max_ttp = Inf; %TODO: Change these values!
handles.pctsettings.min_bbbp = -Inf;
handles.pctsettings.max_bbbp = Inf;
handles.pctsettings.bbbp_first = 1; 
handles.pctsettings.bbbp_last = -1; 
handles.pctsettings.enable_delay = false;

%Display settings
handles.disp.cbv.colormap = 'Jet';
handles.disp.cbv.lo = 0;
handles.disp.cbv.hi = 10;
handles.disp.cbv.scalefactor = 10;

handles.disp.cbf.colormap = 'Jet';
handles.disp.cbf.lo = 0;
handles.disp.cbf.hi = 80;
handles.disp.cbf.scalefactor = 4;

handles.disp.mtt.colormap = 'Jet';
handles.disp.mtt.lo = 0;
handles.disp.mtt.hi = 16;
handles.disp.mtt.scalefactor = 100;

handles.disp.ttp.colormap = 'Jet';
handles.disp.ttp.lo = 0;
handles.disp.ttp.hi = 100;
handles.disp.ttp.scalefactor = 100;

handles.disp.bbbp.colormap = 'Jet';
handles.disp.bbbp.lo = 0;
handles.disp.bbbp.hi = 10;
handles.disp.bbbp.scalefactor = 40;

handles.disp.raw.colormap = 'Gray';
handles.disp.raw.lo = 0;
handles.disp.raw.hi = 100;
handles.disp.raw.scalefactor = 1;

handles.disp.preprocessed.colormap = 'Gray';
handles.disp.preprocessed.lo = 0;
handles.disp.preprocessed.hi = 50;
handles.disp.preprocessed.scalefactor = 1;

handles.disp.residue.colormap = 'Jet';
handles.disp.residue.lo = 0;
handles.disp.residue.hi = 0.01;
handles.disp.residue.scalefactor = 100000;

% Set controls
set(handles.radiobutton_raw,'Enable', 'off');
set(handles.radiobutton_preprocessed,'Enable', 'off');
set(handles.radiobutton_residue,'Enable', 'off');
set(handles.radiobutton_cbv,'Enable', 'off');
set(handles.radiobutton_cbf,'Enable', 'off');
set(handles.radiobutton_mtt,'Enable', 'off');
set(handles.radiobutton_ttp,'Enable', 'off');
set(handles.radiobutton_permeability,'Enable', 'off');
set(handles.radiobutton_patlak,'Enable', 'off');
set(handles.radiobutton_tac,'Enable', 'off');
set(handles.pb_preprocess, 'Enable', 'off');
set(handles.pushbutton_perfusion, 'Enable', 'off');
set(handles.pushbutton_permeability, 'Enable', 'off');

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = Pcttool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_sequence_Callback(hObject, eventdata, handles)
% hObject    handle to slider_sequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if ~handles.haveData
    return
end
%Update edit box and slider
pos = round(get(hObject,'Value'));
set(hObject,'Value',pos);
set(handles.edit_slide,'String',[num2str(pos) ' / ' num2str(handles.nImages)] );
%Update the image
set(handles.image, 'CData', squeeze(handles.imageData(pos,:,:)));
%Update the histogram
update_histogram()
%Update the state
handles.currentImageNo = pos;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider_sequence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_sequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), ...
            get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pb_load.
function pb_load_Callback(hObject, eventdata, handles)
% hObject    handle to pb_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open_mat_file();


function edit_slide_Callback(hObject, eventdata, handles)
% hObject    handle to edit_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function edit_slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
                   get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_loth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_loth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_loth as text
%    str2double(get(hObject,'String')) returns contents of edit_loth as a double
loth = str2double(get(hObject,'String'));
clim = get(handles.axes_image, 'CLim');
hith = clim(2);
set(handles.axes_image, 'CLim', [loth, hith]);


% --- Executes during object creation, after setting all properties.
function edit_loth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_loth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),...
                   get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_hith_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hith (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hith as text
%    str2double(get(hObject,'String')) returns contents of edit_hith as a double
hith = str2double(get(hObject,'String'));
clim = get(handles.axes_image, 'CLim');
loth = clim(1);
set(handles.axes_image, 'CLim', [loth hith]);


% --- Executes during object creation, after setting all properties.
function edit_hith_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_hith (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),...
                   get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_colormap.
function popup_colormap_Callback(hObject, eventdata, handles)
% hObject    handle to popup_colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(hObject,'String');
newmap = contents{get(hObject, 'Value')};
%set(gcf,'Colormap',newmap);
colormap(newmap);


% --- Executes during object creation, after setting all properties.
function popup_colormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),...
                   get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_preprocess.
function pb_preprocess_Callback(hObject, eventdata, handles)
% hObject    handle to pb_preprocess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ticID = tic;
set(gcbf,'Pointer','watch')
pause(0.001); %For some reason this is necessary to update the pointer
guidata(hObject,handles);
preprocess();
computeResidue();
handles = guidata(hObject);
set(gcbf,'Pointer','arrow')
set(handles.pushbutton_perfusion,'Enable','on');
guidata(hObject,handles);
toc(ticID);


function update_histogram()
%Updates the histogram

%Calculate the histogram
handles = guidata(gcbo);
[bins vals] = pct_hist(get(handles.image, 'CData'), 1024);

%Display the histogram
axes(handles.axes_histogram);
stem(bins,vals,'MarkerSize',1);

%Set the axis limits and figure size
rbins = max(bins)-min(bins);
xmin = min(bins) - 0.01*rbins;
xmax = max(bins) + 0.01*rbins;
ymin = 0;
ymax = prctile(vals,98);
if ymax == ymin
    ymax = ymax + 1;
end
axis([xmin xmax ymin ymax]);
set(handles.axes_histogram,'YTick',[]);


function display_tac(x,y)
%Displays a time-attenuation curve for pixel (x,y)
handles = guidata(gcbo);

%Get TAC
tac = pct_tac(handles.imageData,x,y);
axes(handles.axes_tac);

%Display TAC
plot(tac);

%Show position
set(handles.text_tac,'String',['[' num2str(x) ', ' num2str(y) ']']);


function display_patlak(x,y)
%Displays a Patlak plot for pixel (x,y);
handles = guidata(gcbo);
%Get the first and last frame
f = handles.pctsettings.bbbp_first; 
if handles.pctsettings.bbbp_last == -1
    l = handles.nImages;
else
    l = handles.pctsettings.bbbp_last;
end
%Get the patlak plot
if handles.pctsettings.enable_delay
    xx = handles.patlak.xmap(f:l,y,x);
else
    xx = handles.patlak.xmap(f:l);
end

yy = handles.patlak.ymap(f:l,y,x);
%Show the patlak plot
axes(handles.axes_tac);
%plot(xx,yy)
plot(xx,yy,'rs');
%Show the K1 line
slope = handles.patlak.bbbp(y,x) / 60;
tl = slope .* xx;
hold on
plot(xx,tl);
hold off


            
%Show position
set(handles.text_tac,'String',['[' num2str(x) ', ' num2str(y) ']']);
R2 = handles.patlak.r(y,x);
set(handles.text_rsq,'String',num2str(round(100*R2)/100));
guidata(gcbf,handles);


function open_mat_file()
%Open a data set from a .mat file
handles = guidata(gcbf);

[file_name path_name] = uigetfile('*.mat','Load Pcttool image sequence');
if file_name == 0
    return
end
% %Clear old data
% handles.activeData = 'none';
% handles.activeImage = 'none';
% handles.isSequence = false; 
% handles.haveData = false;
% handles.nImages = 0;
% handles.currentImageNo = 1;
% handles.imageData = [];
% handles.mapData = [];
% handles.image_height = 0;
% handles.image_width = 0;
% handles.havePreprocessed = false;
% handles.haveResidueData = false;
% handles.havePerfusion = false;
% handles.cbv.haveData = false;
% handles.cbv.data = [];
% handles.cbf.haveData = false;
% handles.cbf.data = [];
% handles.mtt.haveData = false;
% handles.mtt.data = [];
% handles.ttp.haveData = false;
% handles.ttp.data = [];
% handles.bbbp.haveData = false;
% handles.bbbp.data = [];
% handles.patlak.haveData = false;
% handles.patlak.xmap = [];
% handles.patlak.ymap = [];
% handles.patlak.r = [];
% handles.patlak.bbbp = [];
% %Enable controls
% set(handles.radiobutton_raw,'Enable', 'on');
% set(handles.pb_preprocess,'Enable', 'on');
% set(handles.radiobutton_tac,'Enable', 'on');
handles = reset_ui(handles);
%Get new data
newPatient = load([path_name file_name]);
handles.imageData = newPatient.data;
handles.settings.aif_x = newPatient.aif_x;
handles.settings.aif_y = newPatient.aif_y;
handles.settings.vof_x = newPatient.vof_x;
handles.settings.vof_y = newPatient.vof_y;
handles.settings.dt = newPatient.dt;
if ndims(handles.imageData) ~= 3
    error('Input must be three-dimensional (T x Y x X)');
end
%else: load image into axis
handles.haveData = true;
handles.activeData = 'raw';
handles.isSequence = true;
guidata(gcbf,handles);
dataChanged();


function Tactool_Cursorpos(hObject,eventdata)
% Get the handles
handles = guidata(gcbo);

if ~handles.haveData
    return
end

%Get cursor position
pos = get(handles.axes_image,'CurrentPoint');
x = round(pos(1,1));
y = round(pos(1,2));
if x <= 0 || y <= 0 || x > handles.image_width || y > handles.image_height
    return
end
%Get the pixel value
if handles.isSequence
    val = handles.imageData(handles.currentImageNo,y,x);
else
    val = handles.mapData(y,x);
end
%Update text boxes
set(handles.text_pixel_pos,'String',['[' num2str(x) ', ' num2str(y) ']']);
set(handles.text_value, 'String', num2str(val));


function preprocess()

%Do the preprocessing
handles = guidata(gcbf); 

setActiveData('raw');
%Make sure the raw data is saved to disk
data = handles.imageData;
save('pct_raw.mat', 'data');
clear data

%handles.imageData(handles.imageData == -3024) = 0;

if handles.settings.enable_segment
    [handles.imageData handles.imageMask] = pct_segment(handles.imageData,...
                                    handles.settings.lower_threshold,...
                                    handles.settings.upper_threshold,...
                                    handles.settings.first_frame);
end

if handles.settings.enable_sfilter
    handles.imageData = pct_filter(handles.imageData,...
                                   handles.settings.sfilter_size,...
                                   handles.settings.sfilter_stddev);

end
if handles.settings.enable_tfilter
    handles.imageData = pct_gaussfilter(handles.imageData,...
                                        handles.settings.tfilter_size);
end
if handles.settings.enable_nldf
    handles.imageData = pct_nldif(handles.imageData);
end

if handles.settings.enable_subtractbase
    handles.imageData = pct_subtractbase(handles.imageData,...
                                         handles.settings.first_frame,...
                                         handles.settings.tracer_factor);
end
handles.AIF = pct_tac(handles.imageData, handles.settings.aif_x,...
                      handles.settings.aif_y);
handles.VOF = pct_tac(handles.imageData, handles.settings.vof_x,...
                      handles.settings.vof_y);
if handles.settings.enable_normalizeAIF
    handles.AIF = pct_aifscaling(handles.AIF, handles.VOF);
end
if handles.settings.enable_hematocrit
    handles.imageData = pct_hematocrit(handles.imageData,...
                                       handles.settings.hematocrit);
end

%Apply min and max values (if set)
handles.imageData = pct_truncatevalues(handles.imageData,...
                                       handles.settings.min_preprocessed,...
                                       handles.settings.max_preprocessed);                                      

                                     
handles.havePreprocessed = true;
set(handles.radiobutton_preprocessed,'Enable','on');

preprocessed = handles.imageData;
save 'pct_preprocessed.mat' preprocessed;
clear preprocessed;
handles.activeData = 'preprocessed';
set(handles.radiobutton_preprocessed,'Enable','on');

guidata(gcbf,handles);


function dataChanged()
%To be called whenever new data is displayed in the main axes (axes_image)

handles = guidata(gcbf);
axes(handles.axes_image);

if handles.isSequence
    %Load image
    [handles.nImages handles.image_height handles.image_width] = ...
        size(handles.imageData);
    handles.image = image(squeeze(handles.imageData(1,:,:)), ...
        'CDataMapping', 'scaled');    
    %set slider properties
    set(handles.edit_slide,'Enable','on');
    set(handles.slider_sequence,'Enable','on');
    set(handles.edit_slide,'String',['1 / ' num2str(handles.nImages)]);
    set(handles.slider_sequence,'Min',1);
    set(handles.slider_sequence,'Max',handles.nImages);
    set(handles.slider_sequence,'Value',1);
    step = 1/handles.nImages;
    set(handles.slider_sequence,'SliderStep',[step 10*step]);
else
    %Disable slider and slide text box
    set(handles.edit_slide,'String','');
    set(handles.edit_slide,'Enable','off');
    min = get(handles.slider_sequence,'Min');
    set(handles.slider_sequence,'Value',min);
    set(handles.slider_sequence,'Enable','off');
    %Load image
    handles.image = image(handles.mapData, 'CDataMapping','scaled');
end
handles.currentImageNo = 1;
%Set image to scale
set(handles.axes_image, 'CLimMode', 'manual');
pos = get(handles.axes_image,'Position');
handles.colorbar = colorbar('location', 'EastOutside');
%[h w] = size(get(handles.image,'CData')); %handles.imageData);
set(handles.axes_image,'Position',pos);
%Remove ticks
set(handles.axes_image,'YTick',[]);
set(handles.axes_image,'XTick',[]);
%Establish the TAC callback
set(handles.image, 'ButtonDownFcn', @on_image_click);
guidata(gcbf, handles);
%Set contrast window
setContrast();
%Create histogram
update_histogram();

    
function computeResidue()

handles = guidata(gcbf);

setActiveData('preprocessed');
  
%Truncate the data in time
handles.imageData = pct_truncatetime(handles.imageData,...
                                     handles.settings.first_frame,...
                                     handles.settings.last_frame);
aif = pct_truncatetime(handles.AIF, handles.settings.first_frame,...
                       handles.settings.last_frame);

%Compute the residue functions                   
handles.imageData = cTSVD(aif, handles.imageData,...
                          handles.settings.truncation,...
                          handles.settings.block_m);
%Correct for time-scaling
handles.imageData = (1/handles.settings.dt)*handles.imageData;
                      
%Update the state
handles.activeData = 'residue'; 

set(handles.uipanel_datadisplay, 'SelectedObject', ...
    handles.radiobutton_residue);
set(handles.radiobutton_residue,'Enable','on');
handles.haveResidueData = true;
%Make sure residue is saved to disk
resdata = handles.imageData;
save('pct_residue.mat','resdata');
clear resdata
guidata(gcbf,handles);
dataChanged();


function patlak_permeability()
%Compute the BBBP map, the Patlak plots, and the R^2 map
handles = guidata(gcbf);

if handles.pctsettings.enable_delay
    [bbbp xmap ymap r] = pct_bbbpdc(handles.imageData,...
                                    handles.cbv.data,...
                                    handles.ttp.data,...
                                    handles.AIF,...
                                    handles.settings.dt,...
                                    handles.settings.tissue_density,...
                                    handles.pctsettings.bbbp_first, ...
                                    handles.pctsettings.bbbp_last, ...
                                    handles.imageMask);    
else
    [bbbp xmap ymap r] = pct_bbbp(handles.imageData,...
                                  handles.cbv.data,...
                                  handles.AIF,...
                                  handles.settings.dt,...
                                  handles.settings.tissue_density,...
                                  handles.pctsettings.bbbp_first, ...
                                  handles.pctsettings.bbbp_last, ...
                                  handles.imageMask);    
end
                              



%Apply min and max values (if set)
bbbp = pct_truncatevalues(bbbp,...
    handles.pctsettings.min_bbbp, handles.pctsettings.max_bbbp);
handles.patlak.haveData = true;
handles.patlak.xmap = xmap;
handles.patlak.ymap = ymap;
handles.patlak.r = r;
handles.patlak.bbbp = bbbp;
guidata(gcbf, handles);

function compute_pct()
handles = guidata(gcbf);
if ~handles.haveResidueData
    %something
    return
end
ticID = tic;
set(gcbf,'Pointer','watch');
pause(0.001);
%TODO: Check if map is already computed
if strcmp(handles.activeData,'residue')
    handles.cbv.data = pct_cbv(handles.imageData,...
                               handles.settings.tissue_density,...
                               handles.settings.dt);
    handles.cbf.data = pct_cbf(handles.imageData,...
                               handles.settings.tissue_density,...
                               handles.settings.dt);
    handles.mtt.data = pct_mtt(handles.imageData, handles.settings.dt);
    cmap = importdata('pct_preprocessed.mat');
    handles.ttp.data = pct_ttp(cmap, handles.settings.dt);
else 
    R = importdata('pct_residue.mat');
    handles.cbv.data = pct_cbv(R, handles.settings.tissue_density,...
                               handles.settings.dt);
    handles.cbf.data = pct_cbf(R, handles.settings.tissue_density,...
                               handles.settings.dt);
    handles.mtt.data = pct_mtt(R, handles.settings.dt);    
    if strcmp(handles.activeData,'preprocessed')
        handles.ttp.data = pct_ttp(handles.imageData, handles.settings.dt);
    else
        cmap = importdata('pct_preprocessed.mat');
        handles.ttp.data = pct_ttp(cmap, handles.settings.dt);
    end 
end
%Apply min and max values (if set)
handles.cbv.data = pct_truncatevalues(handles.cbv.data,...
    handles.pctsettings.min_cbv, handles.pctsettings.max_cbv);
handles.cbf.data = pct_truncatevalues(handles.cbf.data,...
    handles.pctsettings.min_cbf, handles.pctsettings.max_cbf);
handles.mtt.data = pct_truncatevalues(handles.mtt.data,...
    handles.pctsettings.min_mtt, handles.pctsettings.max_mtt);
handles.ttp.data = pct_truncatevalues(handles.ttp.data,...
    handles.pctsettings.min_ttp, handles.pctsettings.max_ttp);
%Apply mask
handles.cbv.data = handles.imageMask .* handles.cbv.data;
handles.cbf.data = handles.imageMask .* handles.cbf.data;
handles.mtt.data = handles.imageMask .* handles.mtt.data;
handles.ttp.data = handles.imageMask .* handles.ttp.data;
%Update state
handles.havePerfusion = true;
handles.cbv.haveData = true;
handles.cbf.haveData = true;
handles.mtt.haveData = true;
handles.ttp.haveData = true;
%Update UI controls
set(handles.radiobutton_cbv,'Enable','on');
set(handles.radiobutton_cbf,'Enable','on');
set(handles.radiobutton_mtt,'Enable','on');
set(handles.radiobutton_ttp,'Enable','on');
set(handles.pushbutton_permeability,'Enable','on');
set(gcbf,'Pointer','arrow');
toc(ticID);
guidata(gcbf,handles);


function on_image_click(hObject, eventdata)
%Callback for when a mouse-click is made on the image
%Displays a time-attenuation curve for that pixel

%Get figure handles
handles = guidata(gcbf);
if ~handles.haveData
    return
end

%Get cursor position
pos = get(handles.axes_image,'CurrentPoint');
x = round(pos(1,1));
y = round(pos(1,2));
if x <= 0 || y <= 0 || x > handles.image_width || y > handles.image_height
    return
end

%Save the current (x,y)
handles.x = x;
handles.y = y;
guidata(gcbf,handles);

currentSelection = get(handles.uipanel_plot, 'SelectedObject');
if currentSelection == handles.radiobutton_tac
    display_tac(x,y);
elseif currentSelection == handles.radiobutton_patlak
    display_patlak(x,y);
end

       
function setContrast()
%Sets the contrast window and colormap according to the standard settings
handles = guidata(gcbf);
if handles.isSequence
    switch handles.activeData
        case 'raw'
            ds = handles.disp.raw;
        case 'preprocessed'
            ds = handles.disp.preprocessed;
        case 'residue'
            ds = handles.disp.residue;
    end
else
    switch handles.activeImage
        case 'cbv'
            ds = handles.disp.cbv;
        case 'cbf'
            ds = handles.disp.cbf;
        case 'mtt'
            ds = handles.disp.mtt;
        case 'ttp'
            ds = handles.disp.ttp;
        case 'bbbp'
            ds = handles.disp.bbbp;
    end
end
%Update image
axes(handles.axes_image);
set(handles.axes_image, 'CLim', [ds.lo, ds.hi]);  
colormap(ds.colormap);
%Update ui controls
set(handles.edit_loth,'String',num2str(ds.lo));
set(handles.edit_hith,'String',num2str(ds.hi));
strarray = get(handles.popup_colormap,'String');
ind = strmatch(ds.colormap,strarray);
set(handles.popup_colormap,'Value',ind);
guidata(gcbf,handles);


function [out] = prepareDicom(img)
%prepares an image to be saved as a DICOM file
handles = guidata(gcbf);
if handles.isSequence
    switch handles.activeData
        case 'raw'
            out = int16(img * handles.disp.raw.scalefactor);
        case 'preprocessed'
            out = int16(img * handles.disp.preprocessed.scalefactor);
        case 'residue'
            out = int16(img * handles.disp.raw.scalefactor);
    end
else
    switch handles.activeImage
        case 'cbv'
            out = int16(img * handles.disp.cbv.scalefactor);
        case 'cbf'
            out = int16(img * handles.disp.cbf.scalefactor);
        case 'mtt'
            out = int16(img * handles.disp.mtt.scalefactor);
        case 'ttp'
            out = int16(img * handles.disp.ttp.scalefactor);
        case 'bbbp'
            out = int16(img * handles.disp.bbbp.scalefactor);
    end
end


function setActiveData(newData)
%Sets the active data to newData
%newData: string containing either 'raw','preprocessed', or 'residue'
handles = guidata(gcbf);
inputValid = strcmp(newData,'raw') + strcmp(newData,'preprocessed') + ...
    strcmp(newData,'residue');
if ~inputValid
    disp('Warning: Invalid input to setActiveData');
    return
elseif strcmp(newData,handles.activeData)
    %Just return then
    handles.activeData = newData;
    handles.isSequence = true;
    guidata(gcbf,handles);
    return
end
%Otherwise save the current active data to disk
fileName = ['pct_' handles.activeData '.mat'];
data = handles.imageData;
save(fileName, 'data');
clear data;

newFile = ['pct_' newData '.mat'];
handles.imageData = importdata(newFile);
handles.activeData = newData;
handles.isSequence = true;
guidata(gcbf,handles);

% function open_dicom_directory()
% %Call this function to open a DICOM directory
% handles = guidata(gcbf);
% 
% path_name = uigetdir('.','Load DICOM image sequence');
% if path_name == 0
%     return
% end
% 
% dcm_list = dir([path_name '/*.dcm']);
% if length(dcm_list) < 3
%     %No dicom images
%     return
% end
% %Now that we established that the directory contains DICOM images, proceed to
% %load the header files to determine no. of slices
% dcm_list = dcm_list(3:end);
% %%%FINISH FUNCTION



function [out] = reset_ui(handles)
%Call this function when opening a new file / dataset

%Clear old data
handles.activeData = 'none';
handles.activeImage = 'none';
handles.isSequence = false; 
handles.haveData = false;
handles.nImages = 0;
handles.currentImageNo = 1;
handles.imageData = [];
handles.mapData = [];
handles.image_height = 0;
handles.image_width = 0;
handles.havePreprocessed = false;
handles.haveResidueData = false;
handles.havePerfusion = false;
handles.cbv.haveData = false;
handles.cbv.data = [];
handles.cbf.haveData = false;
handles.cbf.data = [];
handles.mtt.haveData = false;
handles.mtt.data = [];
handles.ttp.haveData = false;
handles.ttp.data = [];
handles.bbbp.haveData = false;
handles.bbbp.data = [];
handles.patlak.haveData = false;
handles.patlak.xmap = [];
handles.patlak.ymap = [];
handles.patlak.r = [];
handles.patlak.bbbp = [];
%Enable controls
set(handles.radiobutton_raw,'Enable', 'on');
set(handles.pb_preprocess,'Enable', 'on');
set(handles.radiobutton_tac,'Enable', 'on');


out = handles;



% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_open_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open_mat_file()


% --------------------------------------------------------------------
function menu_settings_Callback(hObject, eventdata, handles)
% hObject    handle to menu_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_settings_preprocessing_Callback(hObject, eventdata, handles)
% hObject    handle to menu_settings_preprocessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Pctsettings('Pcttool',handles.figure_pcttool);

% --- Executes on button press in pushbutton_test.
function pushbutton_test_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%disp(get(handles.image,'Position'));
assignin('base','handles',handles)
disp('Handles dumped');


% --- Executes when selected object is changed in uipanel_datadisplay.
function uipanel_datadisplay_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_datadisplay 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was 
%             selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
set(gcbf,'Pointer','watch');
pause(0.001);

switch eventdata.NewValue
    case handles.radiobutton_raw
        setActiveData('raw');
        handles = guidata(hObject);
%         fileName = ['pct_' handles.activeData '.mat'];
%         data = handles.imageData;
%         save(fileName, 'data');
%         clear data;
%         handles.imageData = importdata('pct_raw.mat');
%         handles.activeData = 'raw';
%         handles.isSequence = true;
    case handles.radiobutton_preprocessed
        setActiveData('preprocessed');
        handles = guidata(hObject);
%         fileName = ['pct_' handles.activeData '.mat'];
%         data = handles.imageData;
%         save(fileName, 'data');
%         clear data;        
%         handles.imageData = importdata('pct_preprocessed.mat');
%         handles.activeData = 'preprocessed';
%         handles.isSequence = true;
    case handles.radiobutton_residue
        setActiveData('residue');
        handles = guidata(hObject);
%         fileName = ['pct_' handles.activeData '.mat'];
%         data = handles.imageData;
%         save(fileName, 'data');
%         clear data;        
%         handles.imageData = importdata('pct_residue.mat');
%         handles.activeData = 'residue';
%         handles.isSequence = true;
    case handles.radiobutton_cbv
        handles.mapData = handles.cbv.data;
        handles.activeImage = 'cbv';
        handles.isSequence = false;
    case handles.radiobutton_cbf
        handles.mapData = handles.cbf.data;
        handles.activeImage = 'cbf';
        handles.isSequence = false;
    case handles.radiobutton_mtt
        handles.mapData = handles.mtt.data;
        handles.activeImage = 'mtt';
        handles.isSequence = false;
    case handles.radiobutton_ttp
        handles.mapData = handles.ttp.data;
        handles.activeImage = 'ttp';
        handles.isSequence = false;        
    case handles.radiobutton_permeability
        handles.mapData = handles.patlak.bbbp;
        handles.activeImage = 'bbbp'; 
        handles.isSequence = false;
    otherwise
        disp(['Warning: Unknown argument to' ...
            'uipanel_datadisplay_SelectionChangeFcn']);
        
end
guidata(hObject, handles);
dataChanged();
set(gcbf,'Pointer','arrow');


% --- Executes during object creation, after setting all properties.
function uipanel_datadisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_datadisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton_permeability.
function pushbutton_permeability_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_permeability (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.haveData || ~handles.cbv.haveData
    %Notify the user
    return
end
ticID = tic;
set(gcbf,'Pointer','watch');
pause(0.001);

if ~strcmp(handles.activeData,'preprocessed')
    %Get the preprocessed maps
    %handles.activeData = 'preprocessed';
    %set(handles.uipanel_datadisplay, 'SelectedObject', ...
    %handles.radiobutton_preprocessed);
    handles.imageData = importdata('pct_preprocessed.mat');
end

%Calculate the permeability map
guidata(hObject, handles);
patlak_permeability();
handles = guidata(hObject);

%Enable Permeability and Patlak radio buttons
set(handles.radiobutton_permeability,'Enable','on');
set(handles.radiobutton_patlak,'Enable','on');

handles.activeData = 'preprocessed';
set(handles.uipanel_datadisplay, 'SelectedObject', ...
handles.radiobutton_preprocessed);
set(gcbf,'Pointer','arrow');
guidata(hObject,handles);
dataChanged();
toc(ticID);



% --- Executes on button press in pushbutton_perfusion.
function pushbutton_perfusion_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_perfusion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
compute_pct();


% --------------------------------------------------------------------
function menu_settings_perfusion_Callback(hObject, eventdata, handles)
% hObject    handle to menu_settings_perfusion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%  HERE!!!  %%%%%%%%%
Pctflexsettings('Pcttool',handles.figure_pcttool);

% --------------------------------------------------------------------
function menu_settings_permeability_Callback(hObject, eventdata, handles)
% hObject    handle to menu_settings_permeability (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in uipanel_plot.
function uipanel_plot_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_plot 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch eventdata.NewValue
    case handles.radiobutton_tac
        set(handles.text_rsq,'String','0');
        set(handles.text_rsquared,'Visible','off');
        set(handles.text_rsq,'Visible','off');
    case handles.radiobutton_patlak
        set(handles.text_rsquared,'Visible','on');
        set(handles.text_rsq,'Visible','on');        
end


% --------------------------------------------------------------------
function menu_file_save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.haveData
    %Message to user
    return
end

[FileName,PathName,FilterIndex] = uiputfile({'*.dcm'},...
            'Save current image as');
if FileName == 0
    return
end

img = get(handles.image,'CData');
img = prepareDicom(img);
dicomwrite(img,FileName);
        
          
            
% --------------------------------------------------------------------
function menu_settings_display_Callback(hObject, eventdata, handles)
% hObject    handle to menu_settings_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_dumpplot.
function pushbutton_dumpplot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_dumpplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.x == 0 || handles.y == 0
    return
end


currentSelection = get(handles.uipanel_plot, 'SelectedObject');
if currentSelection == handles.radiobutton_tac
    %Dump TAC   
    tac.x = handles.x;
    tac.y = handles.y;
    tac.data = pct_tac(handles.imageData,tac.x,tac.y);
    assignin('base','tac',tac)
elseif currentSelection == handles.radiobutton_patlak
    %Dump Patlak plot
    f = handles.pctsettings.bbbp_first;
    if handles.pctsettings.bbbp_last == -1
        l = handles.nImages;
    else
        l = handles.pctsettings.bbbp_last;
    end
    %Get the patlak plot
    patlak.x = handles.patlak.xmap(f:l,handles.y,handles.x);
    patlak.y = handles.patlak.ymap(f:l,handles.y,handles.x);
    %Get the slope
    patlak.bbbp = handles.patlak.bbbp(handles.y,handles.x) / 60;
    %Get the rsq
    patlak.r2 = handles.patlak.r(handles.y,handles.x);
    assignin('base','patlak',patlak)
end


function varargout = Pctlite(varargin)
% PCTLITE M-file for Pctlite.fig
%      PCTLITE, by itself, creates a new PCTLITE or raises the existing
%      singleton*.
%
%      H = PCTLITE returns the handle to a new PCTLITE or the handle to
%      the existing singleton*.
%
%      PCTLITE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PCTLITE.M with the given input arguments.
%
%      PCTLITE('Property','Value',...) creates a new PCTLITE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Pctlite_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Pctlite_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Pctlite

% Last Modified by GUIDE v2.5 21-Jun-2012 22:13:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Pctlite_OpeningFcn, ...
                   'gui_OutputFcn',  @Pctlite_OutputFcn, ...
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


% --- Executes just before Pctlite is made visible.
function Pctlite_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Pctlite (see VARARGIN)

% Choose default command line output for Pctlite
handles.output = hObject;


%My initialization


axes(handles.axes_image);
set(gca,'DataAspectRatio', [1 1 1]);
set(gcf,'WindowButtonMotionFcn', @Tactool_Cursorpos);

%Image classes
handles.activeData = 'none'; %Could also be 'raw','preprocessed',and 'residue'.
handles.haveData = false;
handles.nImages = 0;
handles.currentImageNo = 1;
handles.imageData = [];

handles.havePreprocessed = false;
handles.haveResidueData = false;

handles.cbv.haveData = false;
handles.cbv.data = 0;
handles.cbf.haveData = false;
handles.cbf.data = 0;
handles.mtt.haveData = false;
handles.mtt.data = 0;
handles.ttp.haveData = false;
handles.ttp.data = 0;
handles.bbbp.haveData = false;
handles.bbbp.data = 0;
handles.patlak.haveData = false;
handles.patlak.data = 0;

%Initialize the settings struct with standard values
handles.settings.tissue_density = 1.05;
handles.settings.hematocrit = 0.73;
handles.settings.lower_threshold = 0;
handles.settings.upper_threshold = 120;
handles.settings.sfilter_size = 5;
handles.settings.sfilter_stddev = 0.5;
handles.settings.tfilter_size = 5;
handles.settings.tfilter_stddev = 0.5;
handles.settings.aif_x = 254; %Change these default values after debugging
handles.settings.aif_y = 206;
handles.settings.vof_x = 237;
handles.settings.vof_y = 411;
handles.settings.truncation = 0.3;
handles.settings.block_m = 2;
handles.settings.first_frame = 2;
handles.settings.last_frame = -1;
handles.settings.tracer_factor = 1;
handles.settings.enable_sfilter = true;
handles.settings.enable_tfilter = true;
handles.settings.enable_segment = true;
handles.settings.enable_subtractbase = true;
handles.settings.enable_normalizeAIF = true;
handles.settings.enable_hematocrit = true;
handles.settings.enable_delay = true;

% Set controls
set(handles.radiobutton_preprocessed,'Enable', 'off');
set(handles.radiobutton_residue,'Enable', 'off');
%set(handles.radiobutton_permeability,'Enable', 'off');

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = Pctlite_OutputFcn(hObject, eventdata, handles) 
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

%Update text boxes
set(handles.text_pixel_pos,'String',['[' num2str(x) ', ' num2str(y) ']']);
set(handles.text_value, 'String',...
    num2str(handles.imageData(handles.currentImageNo,y,x)));


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

display_tac(x,y);


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


function open_mat_file()
%Open a data set from a .mat file

handles = guidata(gcbf);

[file_name path_name] = uigetfile('*.mat','Load image sequence');
handles.imageData = importdata([path_name file_name]);
if ndims(handles.imageData) ~= 3
    error('Input must be three-dimensional (T x Y x X)');
end
%else
%load image into axis
handles.haveData = true;
handles.activeData = 'rawdata';
guidata(gcbf,handles);
dataChanged();


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
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_settings_preprocessing_Callback(hObject, eventdata, handles)
% hObject    handle to menu_settings_preprocessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Pctsettings('Pcttool',handles.figure_pcttool);


function preprocess()

%Do the preprocessing
handles = guidata(gcbf); 

if ~handles.haveData
    return
end

%Save the raw data to disk
rawdata = handles.imageData;
save 'pct_raw.mat' rawdata;
clear rawdata;
%handles.activeData = 'preprocessed';
%delete(handles.image);

if handles.settings.enable_sfilter
    handles.imageData = pct_filter(handles.imageData,...
                                   handles.settings.sfilter_size,...
                                   handles.settings.sfilter_stddev);
end
if handles.settings.enable_tfilter
    handles.imageData = pct_gaussfilter(handles.imageData,...
                                        handles.settings.tfilter_size,...
                                        handles.settings.tfilter_stddev);
end
if handles.settings.enable_segment
    handles.imageData = pct_segment(handles.imageData,...
                                    handles.settings.lower_threshold,...
                                    handles.settings.upper_threshold,...
                                    handles.settings.first_frame);
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
%Truncate the data in time
handles.imageData = pct_truncate(handles.imageData,...
                                 handles.settings.first_frame,...
                                 handles.settings.last_frame);
if handles.settings.last_frame == -1
    handles.AIF = handles.AIF(handles.settings.first_frame:end);
    handles.VOF = handles.VOF(handles.settings.first_frame:end);
else                   
handles.AIF = handles.AIF(handles.settings.first_frame:...
                          handles.settings.last_frame);
handles.VOF = handles.VOF(handles.settings.first_frame:...
                          handles.settings.last_frame);
end
                                     
handles.havePreprocessed = true;
set(handles.radiobutton_preprocessed,'Enable','on');

%Now calculate the residue functions

preprocessed = handles.imageData;
save 'pct_preprocessed.mat' preprocessed;
clear preprocessed;
handles.activeData = 'preprocessed';
set(handles.radiobutton_preprocessed,'Enable','on');
set(handles.uipanel_datadisplay, 'SelectedObject', ...
    handles.radiobutton_preprocessed);

% handles.imageData = cTSVD(handles.AIF,...
%                           handles.imageData,...
%                           handles.settings.truncation,...
%                           handles.settings.block_m);
% 
% handles.activeData = 'residue'; 
% set(handles.radiobutton_residue,'Enable','on');
% set(handles.uipanel_datadisplay, 'SelectedObject', ...
%     handles.radiobutton_residue);

guidata(gcbf,handles);
%dataChanged();



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

fileName = ['pct_' handles.activeData '.mat'];
data = handles.imageData;
save(fileName, 'data');
clear data;

switch eventdata.NewValue
    case handles.radiobutton_raw
        %set(handles.text_test, 'String', 'rdrd');
        handles.imageData = importdata('pct_raw.mat');
        handles.activeData = 'raw';
    case handles.radiobutton_preprocessed
        %set(handles.text_test, 'String', 'ppp');
        handles.imageData = importdata('pct_preprocessed.mat');
        handles.activeData = 'preprocessed';
    case handles.radiobutton_residue
        %set(handles.text_test, 'String', 'rrr');
        handles.imageData = importdata('pct_residue.mat');
        handles.activeData = 'residue';
    otherwise
        %set(handles.text_test, 'String', 'ooo');
end
guidata(hObject, handles);
dataChanged();
set(gcbf,'Pointer','arrow');


% --- Executes during object creation, after setting all properties.
function uipanel_datadisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_datadisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function dataChanged()
%To be called whenever new data is displayed in the main axes (axes_image)

handles = guidata(gcbf);
 
[handles.nImages handles.image_height handles.image_width] = ...
        size(handles.imageData);
set(handles.edit_slide,'String',['1 / ' num2str(handles.nImages)]);

axes(handles.axes_image);
handles.image = image(squeeze(handles.imageData(1,:,:)), ...
        'CDataMapping', 'scaled');
set(handles.axes_image, 'CLimMode', 'manual');
set(handles.axes_image, 'CLim', [0.0, 100.0]);    

handles.colorbar = colorbar('location', 'EastOutside');
%Set image to scale
[t h w] = size(handles.imageData);
pos = get(handles.axes_image,'Position');
set(handles.axes_image,'Position',[pos(1) pos(2) w h]);

%Remove ticks
set(handles.axes_image,'YTick',[]);
set(handles.axes_image,'XTick',[]);

%set slider properties
set(handles.slider_sequence,'Min',1);
set(handles.slider_sequence,'Max',handles.nImages);
set(handles.slider_sequence,'Value',1);
step = 1/handles.nImages;
set(handles.slider_sequence,'SliderStep',[step 10*step]);

%Create histogram
guidata(gcbf, handles);
update_histogram()
handles = guidata(gcbf);

%Establish the TAC callback
%disp(get(handles.axes_image,'HitTest'))
set(handles.image, 'ButtonDownFcn', @on_image_click);
%disp(handles);
guidata(gcbf, handles);
    
function computeResidue()

handles = guidata(gcbf);

handles.imageData = cTSVD(handles.AIF,...
                          handles.imageData,...
                          handles.settings.truncation,...
                          handles.settings.block_m);

handles.activeData = 'residue'; 
set(handles.radiobutton_residue,'Enable','on');
set(handles.uipanel_datadisplay, 'SelectedObject', ...
    handles.radiobutton_residue);

handles.haveResidueData = true;
guidata(gcbf,handles);
dataChanged();
%disp('cheerio!');   



% --- Executes on button press in pushbutton_cbv.
function pushbutton_cbv_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cbv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.haveResidueData
    %something
    return
end
%TODO: Check if map is already computed
if strcmp(handles.activeData,'residue')
    handles.cbv.data = pct_cbv(handles.imageData,...
                              handles.settings.tissue_density);
else
    R = importdata('pct_residue.mat');
    handles.cbv.data = pct_cbv(R, handles.settings.tissue_density);
end
handles.cbv.haveData = true;
guidata(hObject,handles);
imtool(handles.cbv.data,'DisplayRange',[0 10],'Colormap',jet);


% --- Executes on button press in pushbutton_cbf.
function pushbutton_cbf_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cbf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.haveResidueData
    %something
    return
end
%TODO: Check if map is already computed
if strcmp(handles.activeData,'residue')
    handles.cbf.data = pct_cbv(handles.imageData,...
                              handles.settings.tissue_density);
else
    R = importdata('pct_residue.mat');
    handles.cbf.data = pct_cbf(R, handles.settings.tissue_density);
end
handles.cbf.haveData = true;
guidata(hObject,handles);
imtool(handles.cbf.data,'DisplayRange',[0 80],'Colormap',jet);


% --- Executes on button press in pushbutton_mtt.
function pushbutton_mtt_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_mtt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.haveResidueData
    %something
    return
end
%TODO: Check if map is already computed
if strcmp(handles.activeData,'residue')
    handles.mtt.data = pct_mtt(handles.imageData);
else
    R = importdata('pct_residue.mat');
    handles.mtt.data = pct_mtt(R);
end
handles.mtt.haveData = true;
imtool(handles.mtt.data,'DisplayRange',[0 16],'Colormap',jet);


% --- Executes on button press in pushbutton_ttp.
function pushbutton_ttp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ttp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.havePreprocessed
    %something
    return
end
%TODO: Check if map is already computed
if strcmp(handles.activeData,'preprocessed')
    handles.ttp.data = pct_ttp(handles.imageData);
else
    cmap = importdata('pct_preprocessed.mat');
    handles.ttp.data = pct_ttp(cmap);
end
handles.ttp.haveData = true;
guidata(hObject,handles);
imtool(handles.ttp.data,'DisplayRange',[0 handles.nImages],'Colormap',jet);


% % --- Executes on button press in pushbutton10.
% function pushbutton10_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton10 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% 
% % --- Executes on button press in pushbutton_bbbp.
% function pushbutton_bbbp_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton_bbbp (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)



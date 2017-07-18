function varargout = Tactool(varargin)
% TACTOOL M-file for Tactool.fig
%      TACTOOL, by itself, creates a new TACTOOL or raises the existing
%      singleton*.
%
%      H = TACTOOL returns the handle to a new TACTOOL or the handle to
%      the existing singleton*.
%
%      TACTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TACTOOL.M with the given input arguments.
%
%      TACTOOL('Property','Value',...) creates a new TACTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Tactool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Tactool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Tactool

% Last Modified by GUIDE v2.5 18-Jun-2012 12:35:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Tactool_OpeningFcn, ...
                   'gui_OutputFcn',  @Tactool_OutputFcn, ...
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


% --- Executes just before Tactool is made visible.
function Tactool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Tactool (see VARARGIN)

% Choose default command line output for Tactool
handles.output = hObject;

%My initialization
handles.image_loaded = false;
handles.image_no = 1;
axes(handles.axes_image);
set(gca,'DataAspectRatio', [1 1 1]);
set(gcf,'WindowButtonMotionFcn', @Tactool_Cursorpos);

%Take the input argument
if nargin > 0
    handles.image_data = varargin{1};
    %Temporary solution
    if ndims(handles.image_data) ~= 3
        error('Input must be three-dimensional (T x Y x X)');
    end
    %else
    %load image into axis
    handles.image_loaded = true;
    [handles.image_length handles.image_height handles.image_width] = size(handles.image_data);
    set(handles.edit_slide,'String',['1 / ' num2str(handles.image_length)]);
    axes(handles.axes_image);
    handles.image = image(squeeze(handles.image_data(1,:,:)), 'CDataMapping', 'scaled');
    set(handles.axes_image, 'CLimMode', 'manual');
    set(handles.axes_image, 'CLim', [0.0, 100.0]);
    
    handles.colorbar = colorbar('location', 'EastOutside');
    %Set image to scale
    [t h w] = size(handles.image_data);
    pos = get(handles.axes_image,'Position');
    set(handles.axes_image,'Position',[pos(1) pos(2) w h]);
    
    %Remove ticks
    set(handles.axes_image,'YTick',[]);
    set(handles.axes_image,'XTick',[]);
    
    %set slider properties
    set(handles.slider_sequence,'Min',1);
    set(handles.slider_sequence,'Max',handles.image_length);
    set(handles.slider_sequence,'Value',1);
    step = 1/handles.image_length;
    set(handles.slider_sequence,'SliderStep',[step 10*step]);
    
    %Create histogram
    guidata(hObject, handles);
    update_histogram()
    
    %Establish the TAC callback
    disp(get(handles.axes_image,'HitTest'))
    set(handles.image, 'ButtonDownFcn', @on_image_click);
    
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Tactool wait for user response (see UIRESUME)
% uiwait(handles.figure_tactool);


% --- Outputs from this function are returned to the command line.
function varargout = Tactool_OutputFcn(hObject, eventdata, handles) 
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
if ~handles.image_loaded
    return
end
%Update edit box and slider
pos = round(get(hObject,'Value'));
set(hObject,'Value',pos);
set(handles.edit_slide,'String',[num2str(pos) ' / ' num2str(handles.image_length)] );
%Update the image
set(handles.image, 'CData', squeeze(handles.image_data(pos,:,:)));
%Update the histogram
update_histogram()


% --- Executes during object creation, after setting all properties.
function slider_sequence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_sequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in pb_load.
function pb_load_Callback(hObject, eventdata, handles)
% hObject    handle to pb_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file_name path_name] = uigetfile('*.mat','Load image sequence');
handles.image_data = importdata([path_name file_name]);
if ndims(handles.image_data) ~= 3
    error('Input must be three-dimensional (T x Y x X)');
end
%else
%load image into axis
handles.image_loaded = true;
[handles.image_length handles.image_height handles.image_width] = size(handles.image_data);
set(handles.edit_slide,'String',['1 / ' num2str(handles.image_length)]);
axes(handles.axes_image);
handles.image = image(squeeze(handles.image_data(1,:,:)), 'CDataMapping', 'scaled');
set(handles.axes_image, 'CLimMode', 'manual');
set(handles.axes_image, 'CLim', [0.0, 100.0]);

handles.colorbar = colorbar('location', 'EastOutside');
%Set image to scale
[t h w] = size(handles.image_data);
pos = get(handles.axes_image,'Position');
set(handles.axes_image,'Position',[pos(1) pos(2) w h]);

%Remove ticks
set(handles.axes_image,'YTick',[]);
set(handles.axes_image,'XTick',[]);

%set slider properties
set(handles.slider_sequence,'Min',1);
set(handles.slider_sequence,'Max',handles.image_length);
set(handles.slider_sequence,'Value',1);
step = 1/handles.image_length;
set(handles.slider_sequence,'SliderStep',[step 10*step]);

%Create histogram
guidata(hObject, handles);
update_histogram()

%Establish the TAC callback
disp(get(handles.axes_image,'HitTest'))
set(handles.image, 'ButtonDownFcn', @on_image_click);




function edit_slide_Callback(hObject, eventdata, handles)
% hObject    handle to edit_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_slide as text
%        str2double(get(hObject,'String')) returns contents of edit_slide as a double


% --- Executes during object creation, after setting all properties.
function edit_slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Tactool_Cursorpos(hObject,eventdata)
% Get the handles
handles = guidata(gcbo);

if ~handles.image_loaded
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
set(handles.text_value, 'String',num2str(handles.image_data(handles.image_no,x,y)));



function edit_loth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_loth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_loth as text
%        str2double(get(hObject,'String')) returns contents of edit_loth as a double
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
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_hith_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hith (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hith as text
%        str2double(get(hObject,'String')) returns contents of edit_hith as a double
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
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_colormap.
function popup_colormap_Callback(hObject, eventdata, handles)
% hObject    handle to popup_colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popup_colormap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_colormap
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
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_histogram.
function pb_histogram_Callback(hObject, eventdata, handles)
% hObject    handle to pb_histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Compute histogram
[bins vals] = pct_hist(get(handles.image, 'CData'), 1024);

%Display histogram in new window
%handles.histogram = figure('Name','Histogram');
axes(handles.axes_histogram);
stem(bins,vals,'MarkerSize',1);

%Set the axis limits and figure size
rbins = max(bins)-min(bins);
xmin = min(bins) - 0.01*rbins;
xmax = max(bins) + 0.01*rbins;
ymin = 0;
ymax = prctile(vals,98);
axis([xmin xmax ymin ymax]);
set(handles.axes_histogram,'YTick',[]);
%set(handles.histogram,'Units','pixels');
%pos = get(handles.histogram,'Position');
%set(handles.histogram,'Position',[pos(1) pos(2) 690 130]);
%set(handles.histogram,'CloseRequestFcn', @close_histogram);
%set(handles.histogram,'Parent',gcbo);

guidata(hObject, handles);

% function close_histogram(hObject, eventdata)
% %Callback for when the histogram figure is closed
% 
% %histofigure = guidata(gcbo);
% handles = get(gcbo,'Parent');
% disp(get(handles,'Name'));


function update_histogram()
%Updates the histogram

%Calculate the histogram
handles = guidata(gcf);
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
axis([xmin xmax ymin ymax]);
set(handles.axes_histogram,'YTick',[]);


function on_image_click(hObject, eventdata)
%Callback for when a mouse-click is made on the image
%Displays a time-attenuation curve for that pixel

%Get figure handles
handles = guidata(gcbo);
if ~handles.image_loaded
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
tac = pct_tac(handles.image_data,x,y);
axes(handles.axes_tac);

%Display TAC
plot(tac);

%Show position
set(handles.text_tac,'String',['[' num2str(x) ', ' num2str(y) ']']);


%Features to add:
%   -lock histogram axes
%   -superimpose TACs
%   -Zoom/Pan
%   -Set TAC/histogram aces






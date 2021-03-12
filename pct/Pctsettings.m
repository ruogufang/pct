function varargout = pctsettings(varargin)
% PCTSETTINGS M-file for pctsettings.fig
%      PCTSETTINGS, by itself, creates a new PCTSETTINGS or raises the existing
%      singleton*.
%
%      H = PCTSETTINGS returns the handle to a new PCTSETTINGS or the handle to
%      the existing singleton*.
%
%      PCTSETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PCTSETTINGS.M with the given input arguments.
%
%      PCTSETTINGS('Property','Value',...) creates a new PCTSETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pctsettings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pctsettings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pctsettings

% Last Modified by GUIDE v2.5 02-Jul-2012 13:45:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pctsettings_OpeningFcn, ...
                   'gui_OutputFcn',  @pctsettings_OutputFcn, ...
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


% --- Executes just before pctsettings is made visible.
function pctsettings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pctsettings (see VARARGIN)

% Choose default command line output for pctsettings
handles.output = hObject;

%Check to see if Tacttool was passed as parameter
mainGuiInput = find(strcmp(varargin, 'Pcttool'));
if (isempty(mainGuiInput)) || (length(varargin) <= mainGuiInput) || (~ishandle(varargin{mainGuiInput+1}))
    handles.dontOpen = true;
else
    % Remember the handle
    handles.mainGUI = varargin{mainGuiInput+1};
    mainGuiHandles = guidata(handles.mainGUI);
    handles.settings = mainGuiHandles.settings;
end

% Update handles structure
guidata(hObject, handles);

%Load the values into the ui elements
load_values(handles.figure_pctsettings);



% --- Outputs from this function are returned to the command line.
function varargout = pctsettings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = [];



% --- Executes when user attempts to close figure.
function pctsettings_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%uiresume(hObject);


% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Save the settings
mainHandles = guidata(handles.mainGUI);
mainHandles.settings = handles.settings;
guidata(handles.mainGUI,mainHandles);

%disp(handles);
%disp(handles.settings);

%close the figure
close(gcbf);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcbf);


function edit_aif_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_aif_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_aif_x as text
%        str2double(get(hObject,'String')) returns contents of edit_aif_x as a double
handles.settings.aif_x = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_aif_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_aif_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_vof_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vof_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vof_x as text
%        str2double(get(hObject,'String')) returns contents of edit_vof_x as a double
handles.settings.vof_x = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_vof_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vof_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_aif_y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_aif_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_aif_y as text
%        str2double(get(hObject,'String')) returns contents of edit_aif_y as a double
handles.settings.aif_y = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_aif_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_aif_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vof_y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vof_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vof_y as text
%        str2double(get(hObject,'String')) returns contents of edit_vof_y as a double
handles.settings.vof_y = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_vof_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vof_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox7.
%function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7



% --- Executes on button press in checkbox_sfilter.
function checkbox_sfilter_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_sfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_sfilter
handles.settings.enable_sfilter = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in checkbox_tfilter.
function checkbox_tfilter_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_tfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_tfilter
handles.settings.enable_tfilter = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in checkbox_segment.
function checkbox_segment_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_segment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_segment
handles.settings.enable_segment = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in checkbox_subtractbase.
function checkbox_subtractbase_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_subtractbase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_subtractbase
handles.settings.enable_subtractbase = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in checkbox_normalizeAIF.
function checkbox_normalizeAIF_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_normalizeAIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_normalizeAIF
handles.settings.enable_normalizeAIF = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in checkbox_hematocrit.
function checkbox_hematocrit_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_hematocrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_hematocrit
handles.settings.enable_hematocrit = get(hObject,'Value');
guidata(hObject,handles);


function edit_tissue_density_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tissue_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tissue_density as text
%        str2double(get(hObject,'String')) returns contents of edit_tissue_density as a double
handles.settings.tissue_density = str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_tissue_density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tissue_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_hematocrit_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hematocrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hematocrit as text
%        str2double(get(hObject,'String')) returns contents of edit_hematocrit as a double
handles.settings.hematocrit = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_hematocrit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_hematocrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tracer_factor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tracer_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tracer_factor as text
%        str2double(get(hObject,'String')) returns contents of edit_tracer_factor as a double
handles.settings.tracer_factor = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_tracer_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tracer_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_block_m_Callback(hObject, eventdata, handles)
% hObject    handle to edit_block_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_block_m as text
%        str2double(get(hObject,'String')) returns contents of edit_block_m as a double
handles.settings.block_m = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_block_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_block_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_truncation_Callback(hObject, eventdata, handles)
% hObject    handle to edit_truncation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_truncation as text
%        str2double(get(hObject,'String')) returns contents of edit_truncation as a double
handles.settings.truncation = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_truncation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_truncation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tfilter_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tfilter_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tfilter_size as text
%        str2double(get(hObject,'String')) returns contents of edit_tfilter_size as a double
handles.settings.tfilter_size = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_tfilter_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tfilter_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tfilter_stddev_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tfilter_stddev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tfilter_stddev as text
%        str2double(get(hObject,'String')) returns contents of edit_tfilter_stddev as a double
handles.settings.tfilter_stddev = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_tfilter_stddev_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tfilter_stddev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sfilter_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sfilter_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sfilter_size as text
%        str2double(get(hObject,'String')) returns contents of edit_sfilter_size as a double
handles.settings.sfilter_size = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_sfilter_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sfilter_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sfilter_stddev_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sfilter_stddev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sfilter_stddev as text
%        str2double(get(hObject,'String')) returns contents of edit_sfilter_stddev as a double
handles.settings.sfilter_stddev = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_sfilter_stddev_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sfilter_stddev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_first_frame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_first_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_first_frame as text
%        str2double(get(hObject,'String')) returns contents of edit_first_frame as a double
handles.settings.first_frame = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_first_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_first_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_last_frame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_last_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_last_frame as text
%        str2double(get(hObject,'String')) returns contents of edit_last_frame as a double
handles.settings.last_frame = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_last_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_last_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_lower_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lower_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lower_threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_lower_threshold as a double
handles.settings.lower_threshold = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_lower_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lower_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_upper_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_upper_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_upper_threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_upper_threshold as a double
handles.settings.upper_threshold = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_upper_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_upper_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function load_values(current_figure)
%Loads the variables from the struct 'settings' into the ui elements.
handles = guidata(current_figure);

set(handles.edit_tissue_density, 'String', num2str(handles.settings.tissue_density));
set(handles.edit_hematocrit, 'String', num2str(handles.settings.hematocrit));
set(handles.edit_lower_threshold, 'String', num2str(handles.settings.lower_threshold));
set(handles.edit_upper_threshold, 'String', num2str(handles.settings.upper_threshold));
set(handles.edit_sfilter_size, 'String', num2str(handles.settings.sfilter_size));
set(handles.edit_sfilter_stddev, 'String', num2str(handles.settings.sfilter_stddev));
set(handles.edit_tfilter_size, 'String', num2str(handles.settings.tfilter_size));
set(handles.edit_tfilter_stddev, 'String', num2str(handles.settings.tfilter_stddev));
set(handles.edit_aif_x, 'String', num2str(handles.settings.aif_x));
set(handles.edit_aif_y, 'String', num2str(handles.settings.aif_y));
set(handles.edit_vof_x, 'String', num2str(handles.settings.vof_x));
set(handles.edit_vof_y, 'String', num2str(handles.settings.vof_y));
set(handles.edit_truncation, 'String', num2str(handles.settings.truncation));
set(handles.edit_block_m, 'String', num2str(handles.settings.block_m));
set(handles.edit_first_frame, 'String', num2str(handles.settings.first_frame));
set(handles.edit_last_frame, 'String', num2str(handles.settings.last_frame));
set(handles.edit_tracer_factor, 'String', num2str(handles.settings.tracer_factor));
set(handles.checkbox_sfilter, 'Value', handles.settings.enable_sfilter);
set(handles.checkbox_tfilter, 'Value', handles.settings.enable_tfilter);
set(handles.checkbox_segment, 'Value', handles.settings.enable_segment);
set(handles.checkbox_subtractbase, 'Value', handles.settings.enable_subtractbase);
set(handles.checkbox_normalizeAIF, 'Value', handles.settings.enable_normalizeAIF);
set(handles.checkbox_hematocrit, 'Value', handles.settings.enable_hematocrit);
set(handles.checkbox_nldf, 'Value', handles.settings.enable_nldf);
%set(handles.checkbox_delay, 'Value', handles.settings.enable_delay);


% --- Executes on button press in checkbox_nldf.
function checkbox_nldf_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_nldf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.settings.enable_nldf = get(hObject,'Value');
guidata(hObject,handles);



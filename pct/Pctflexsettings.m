function varargout = Pctflexsettings(varargin)
% PCTFLEXSETTINGS M-file for Pctflexsettings.fig
%      PCTFLEXSETTINGS, by itself, creates a new PCTFLEXSETTINGS or raises the existing
%      singleton*.
%
%      H = PCTFLEXSETTINGS returns the handle to a new PCTFLEXSETTINGS or the handle to
%      the existing singleton*.
%
%      PCTFLEXSETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PCTFLEXSETTINGS.M with the given input arguments.
%
%      PCTFLEXSETTINGS('Property','Value',...) creates a new PCTFLEXSETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Pctflexsettings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Pctflexsettings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Pctflexsettings

% Last Modified by GUIDE v2.5 02-Jul-2012 14:18:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Pctflexsettings_OpeningFcn, ...
                   'gui_OutputFcn',  @Pctflexsettings_OutputFcn, ...
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


% --- Executes just before Pctflexsettings is made visible.
function Pctflexsettings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Pctflexsettings (see VARARGIN)

% Choose default command line output for Pctflexsettings
handles.output = hObject;

%Check to see if Tacttool was passed as parameter
mainGuiInput = find(strcmp(varargin, 'Pcttool'));
if (isempty(mainGuiInput)) || (length(varargin) <= mainGuiInput) || (~ishandle(varargin{mainGuiInput+1}))
    handles.dontOpen = true;
else
    % Remember the handle
    handles.mainGUI = varargin{mainGuiInput+1};
    mainGuiHandles = guidata(handles.mainGUI);
    handles.settings = mainGuiHandles.pctsettings;
end

% Update handles structure
guidata(hObject, handles);






% --- Outputs from this function are returned to the command line.
function varargout = Pctflexsettings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox_main.
function listbox_main_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_main contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_main
%disp(get(hObject,'Value'));
sel = get(hObject,'Value');
load_setting(sel);


% --- Executes during object creation, after setting all properties.
function listbox_main_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_main_Callback(hObject, eventdata, handles)
% hObject    handle to edit_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_main as text
%        str2double(get(hObject,'String')) returns contents of edit_main as a double
val = str2double(get(hObject,'String'));
save_setting(val);


% --- Executes during object creation, after setting all properties.
function edit_main_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Save the settings
mainHandles = guidata(handles.mainGUI);
mainHandles.pctsettings = handles.settings;
guidata(handles.mainGUI,mainHandles);
%close the figure
close(gcbf);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcbf);

function load_setting(num)
handles = guidata(gcbf);
switch num
    case 1
        set(handles.edit_main,'String',num2str(handles.settings.min_cbv));
    case 2
        set(handles.edit_main,'String',num2str(handles.settings.max_cbv));
    case 3
        set(handles.edit_main,'String',num2str(handles.settings.min_cbf));
    case 4
        set(handles.edit_main,'String',num2str(handles.settings.max_cbf));
    case 5
        set(handles.edit_main,'String',num2str(handles.settings.min_mtt));
    case 6
        set(handles.edit_main,'String',num2str(handles.settings.max_mtt));
    case 7
        set(handles.edit_main,'String',num2str(handles.settings.min_ttp));
    case 8
        set(handles.edit_main,'String',num2str(handles.settings.max_ttp));
    case 9
        set(handles.edit_main,'String',num2str(handles.settings.min_bbbp));
    case 10
        set(handles.edit_main,'String',num2str(handles.settings.max_bbbp));
    case 11
        set(handles.edit_main,'String',num2str(handles.settings.bbbp_first));
    case 12
        set(handles.edit_main,'String',num2str(handles.settings.bbbp_last));
end
guidata(gcbf,handles);


function save_setting(val)
handles = guidata(gcbf);
num = get(handles.listbox_main,'Value');
switch num
    case 1
        handles.settings.min_cbv = val;
    case 2
        handles.settings.max_cbv = val;
    case 3
        handles.settings.min_cbf = val;
    case 4
        handles.settings.max_cbf = val;
    case 5
        handles.settings.min_mtt = val;
    case 6
        handles.settings.max_mtt = val;
    case 7
        handles.settings.min_ttp = val;
    case 8
        handles.settings.max_ttp = val;
    case 9
        handles.settings.min_bbbp = val;
    case 10
        handles.settings.max_bbbp = val;
    case 11
        handles.settings.bbbp_first = val;
    case 12
        handles.settings.bbbp_last = val;

end
guidata(gcbf,handles);












function varargout = V4_GUI(varargin)
% V4_GUI MATLAB code for V4_GUI.fig
%      V4_GUI, by itself, creates a new V4_GUI or raises the existing
%      singleton*.
%
%      H = V4_GUI returns the handle to a new V4_GUI or the handle to
%      the existing singleton*.
%
%      V4_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in V4_GUI.M with the given input arguments.
%
%      V4_GUI('Property','Value',...) creates a new V4_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before V4_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to V4_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help V4_GUI

% Last Modified by GUIDE v2.5 21-Mar-2014 19:10:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @V4_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @V4_GUI_OutputFcn, ...
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


% --- Executes just before V4_GUI is made visible.
function V4_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to V4_GUI (see VARARGIN)

% Choose default command line output for V4_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes V4_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

    GUI_tables();

    set(handles.field_table,'Data',field_table);
    set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
    set(handles.ducts_size_table,'Data',ducts_size_table);
    set(handles.insulation_table,'Data',insulation_table);
    set(handles.HEX_table,'Data',HEX_table);
    set(handles.plant_table,'Data',plant_table);
    % optimization:
    set(handles.field_Xb_table,'Data',field_Xb_table);
    set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
    set(handles.constraints_table,'Data',constraints_table);
    
    set(handles.model_type,'Value',model_type+1);
    set(handles.location,'Value',location+1);
    set(handles.plant_type,'Value',plant_type+1);
    set(handles.HTS_type,'Value',HTS_type+1);
    set(handles.system_pressure,'Value',system_pressure+1);
    set(handles.file_name,'String','');




% --- Outputs from this function are returned to the command line.
function varargout = V4_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in location.
function location_Callback(hObject, eventdata, handles)
% hObject    handle to location (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns location contents as cell array
%        contents{get(hObject,'Value')} returns selected item from location
global location
location = get(handles.location,'Value')-1

% --- Executes during object creation, after setting all properties.
function location_CreateFcn(hObject, eventdata, handles)
% hObject    handle to location (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plant_type.
function plant_type_Callback(hObject, eventdata, handles)
% hObject    handle to plant_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plant_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plant_type
global plant_type
plant_type = get(handles.plant_type,'Value')-1

% --- Executes during object creation, after setting all properties.
function plant_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plant_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in HTS_type.
function HTS_type_Callback(hObject, eventdata, handles)
% hObject    handle to HTS_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns HTS_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from HTS_type
global system_pressure...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

% --- Executes during object creation, after setting all properties.
function HTS_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HTS_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in system_pressure.
function system_pressure_Callback(hObject, eventdata, handles)
% hObject    handle to system_pressure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns system_pressure contents as cell array
%        contents{get(hObject,'Value')} returns selected item from system_pressure
global system_pressure
system_pressure = get(handles.system_pressure,'Value')-1

% --- Executes during object creation, after setting all properties.
function system_pressure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to system_pressure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in model_type.
function model_type_Callback(hObject, eventdata, handles)
% hObject    handle to model_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns model_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from model_type
global model_type
model_type = get(handles.model_type,'Value')-1



% --- Executes during object creation, after setting all properties.
function model_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in constraints_table.
function constraints_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to constraints_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected cell(s) is changed in constraints_table.
function constraints_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to constraints_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function constraints_table_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to constraints_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on constraints_table and none of its controls.
function constraints_table_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to constraints_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when entered data in editable cell(s) in insulation_table.
function insulation_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to insulation_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function HTS_Xb_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HTS_Xb_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function dishNreceiver_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dishNreceiver_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in clear_ducts.
function clear_ducts_Callback(hObject, eventdata, handles)
% hObject    handle to clear_ducts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    GUI_tables();   

%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
    set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization:
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in default_ducts.
function default_ducts_Callback(hObject, eventdata, handles)
% hObject    handle to default_ducts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    default();

%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
    set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization V
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in clear_macro.
function clear_macro_Callback(hObject, eventdata, handles)
% hObject    handle to clear_macro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    GUI_tables(); 

%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization:
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
    
    set(handles.model_type,'Value',model_type+1);
    set(handles.location,'Value',location+1);
    set(handles.plant_type,'Value',plant_type+1);
    set(handles.HTS_type,'Value',HTS_type+1);
    set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in default_macro.
function default_macro_Callback(hObject, eventdata, handles)
% hObject    handle to default_macro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    default();
 
%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization V
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
    set(handles.model_type,'Value',model_type+1);
    set(handles.location,'Value',location+1);
    set(handles.plant_type,'Value',plant_type+1);
    set(handles.HTS_type,'Value',HTS_type+1);
    set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in clear_insulation.
function clear_insulation_Callback(hObject, eventdata, handles)
% hObject    handle to clear_insulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    GUI_tables();
        
%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
    set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization:
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in default_insulation.
function default_insulation_Callback(hObject, eventdata, handles)
% hObject    handle to default_insulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    default();
        
%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
    set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization V
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in clear_field.
function clear_field_Callback(hObject, eventdata, handles)
% hObject    handle to clear_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    GUI_tables();

    set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization:
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in default_field.
function default_field_Callback(hObject, eventdata, handles)
% hObject    handle to default_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    default();
  
    set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization V
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in clear_dishNreceiver.
function clear_dishNreceiver_Callback(hObject, eventdata, handles)
% hObject    handle to clear_dishNreceiver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    GUI_tables();

    global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    HTS_type...
    system_pressure...
    file_name
    

%     set(handles.field_table,'Data',field_table);
    set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization:
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in default_dishNreceiver.
function default_dishNreceiver_Callback(hObject, eventdata, handles)
% hObject    handle to default_dishNreceiver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    default();

    global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type
  
%     set(handles.field_table,'Data',field_table);
    set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization V
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in clear_plant.
function clear_plant_Callback(hObject, eventdata, handles)
% hObject    handle to clear_plant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    GUI_tables();

    global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    HTS_type...
    system_pressure...
    file_name
    

%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
    set(handles.plant_table,'Data',plant_table);
%     % optimization:
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in default_plant.
function default_plant_Callback(hObject, eventdata, handles)
% hObject    handle to default_plant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    default();

    global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type
  
%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
    set(handles.plant_table,'Data',plant_table);
%     % optimization V
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in clear_HEX.
function clear_HEX_Callback(hObject, eventdata, handles)
% hObject    handle to clear_HEX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    GUI_tables();

    global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    HTS_type...
    system_pressure...
    file_name
    

%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
    set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization:
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in default_HEX.
function default_HEX_Callback(hObject, eventdata, handles)
% hObject    handle to default_HEX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    default();

    global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type
  
%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
    set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization V
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in clear_project.
function clear_project_Callback(hObject, eventdata, handles)
% hObject    handle to clear_project (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

GUI_tables();
    
    set(handles.field_table,'Data',field_table);
    set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
    set(handles.ducts_size_table,'Data',ducts_size_table);
    set(handles.insulation_table,'Data',insulation_table);
    set(handles.HEX_table,'Data',HEX_table);
    set(handles.plant_table,'Data',plant_table);
%     % optimization:
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
    
    set(handles.model_type,'Value',model_type+1);
    set(handles.location,'Value',location+1);
    set(handles.plant_type,'Value',plant_type+1);
    set(handles.HTS_type,'Value',HTS_type+1);
    set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in default_project.
function default_project_Callback(hObject, eventdata, handles)
% hObject    handle to default_project (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    default();
  
    set(handles.field_table,'Data',field_table);
    set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
    set(handles.ducts_size_table,'Data',ducts_size_table);
    set(handles.insulation_table,'Data',insulation_table);
    set(handles.HEX_table,'Data',HEX_table);
    set(handles.plant_table,'Data',plant_table);
%     % optimization V
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
    
    set(handles.model_type,'Value',model_type+1);
    set(handles.location,'Value',location+1);
    set(handles.plant_type,'Value',plant_type+1);
    set(handles.HTS_type,'Value',HTS_type+1);
    set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in clear_constraints.
function clear_constraints_Callback(hObject, eventdata, handles)
% hObject    handle to clear_constraints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    GUI_tables();
    
%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization:
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
    set(handles.constraints_table,'Data',constraints_table);
    
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in default_constraints.
function default_constraints_Callback(hObject, eventdata, handles)
% hObject    handle to default_constraints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    default();

%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization V
%     set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
    set(handles.constraints_table,'Data',constraints_table);
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in clear_field_Xb.
function clear_field_Xb_Callback(hObject, eventdata, handles)
% hObject    handle to clear_field_Xb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    GUI_tables();

%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization:
    set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in defauult_field_Xb.
function defauult_field_Xb_Callback(hObject, eventdata, handles)
% hObject    handle to defauult_field_Xb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    default();
  
%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization V
    set(handles.field_Xb_table,'Data',field_Xb_table);
%     set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in clear_HTS_Xb.
function clear_HTS_Xb_Callback(hObject, eventdata, handles)
% hObject    handle to clear_HTS_Xb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    GUI_tables();

%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization:
%     set(handles.field_Xb_table,'Data',field_Xb_table);
    set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in default_HTS_Xb.
function default_HTS_Xb_Callback(hObject, eventdata, handles)
% hObject    handle to default_HTS_Xb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    default();

%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
    % optimization V
%     set(handles.field_Xb_table,'Data',field_Xb_table);
    set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
%     set(handles.constraints_table,'Data',constraints_table);
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in clear_optimization.
function clear_optimization_Callback(hObject, eventdata, handles)
% hObject    handle to clear_optimization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    GUI_tables();

%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
%     % optimization:
    set(handles.field_Xb_table,'Data',field_Xb_table);
    set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
    set(handles.constraints_table,'Data',constraints_table);
    
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in default_optimization.
function default_optimization_Callback(hObject, eventdata, handles)
% hObject    handle to default_optimization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    default();

%     set(handles.field_table,'Data',field_table);
%     set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
%     set(handles.ducts_size_table,'Data',ducts_size_table);
%     set(handles.insulation_table,'Data',insulation_table);
%     set(handles.HEX_table,'Data',HEX_table);
%     set(handles.plant_table,'Data',plant_table);
    % optimization V
    set(handles.field_Xb_table,'Data',field_Xb_table);
    set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
    set(handles.constraints_table,'Data',constraints_table);
%     set(handles.model_type,'Value',model_type+1);
%     set(handles.location,'Value',location+1);
%     set(handles.plant_type,'Value',plant_type+1);
%     set(handles.HTS_type,'Value',HTS_type+1);
%     set(handles.system_pressure,'Value',system_pressure+1);
%     set(handles.file_name,'String',file_name);

% --- Executes on button press in clear_all.
function clear_all_Callback(hObject, eventdata, handles)
% hObject    handle to clear_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    HTS_type...
    file_name

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    GUI_tables();

    set(handles.field_table,'Data',field_table);
    set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
    set(handles.ducts_size_table,'Data',ducts_size_table);
    set(handles.insulation_table,'Data',insulation_table);
    set(handles.HEX_table,'Data',HEX_table);
    set(handles.plant_table,'Data',plant_table);
    % optimization:
    set(handles.field_Xb_table,'Data',field_Xb_table);
    set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
    set(handles.constraints_table,'Data',constraints_table);
    
    set(handles.model_type,'Value',model_type+1);
    set(handles.location,'Value',location+1);
    set(handles.plant_type,'Value',plant_type+1);
    set(handles.HTS_type,'Value',HTS_type+1);
    set(handles.system_pressure,'Value',system_pressure+1);
    set(handles.file_name,'String',file_name);

% --- Executes on button press in default_all.
function default_all_Callback(hObject, eventdata, handles)
% hObject    handle to default_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    file_name...
    HTS_type

system_pressure = get(handles.system_pressure,'Value')-1;
HTS_type = get(handles.HTS_type,'Value')-1;

    default();
  
    set(handles.field_table,'Data',field_table);
    set(handles.dishNreceiver_table,'Data',dishNreceiver_table);
    set(handles.ducts_size_table,'Data',ducts_size_table);
    set(handles.insulation_table,'Data',insulation_table);
    set(handles.HEX_table,'Data',HEX_table);
    set(handles.plant_table,'Data',plant_table);
    % optimization V
    set(handles.field_Xb_table,'Data',field_Xb_table);
    set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
    set(handles.constraints_table,'Data',constraints_table);
    
    set(handles.model_type,'Value',model_type+1);
    set(handles.location,'Value',location+1);
    set(handles.plant_type,'Value',plant_type+1);
    set(handles.HTS_type,'Value',HTS_type+1);
    set(handles.system_pressure,'Value',system_pressure+1);
    set(handles.file_name,'String',file_name);

% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global  INPUT...
    field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    HTS_type...
    file_name



p = uigetfile('*.mat','Load Input Data')
k = importdata(p)
INPUT = k.INPUT


        %% Macro Variables:
        
        set(handles.site_location,'Value',INPUT.site_location+1);
        set(handles.plant_type,'Value',INPUT.plant_type+1);
        set(handles.system_pressure,'Value',INPUT.system_pressure+1);
        set(handles.model_type,'Value',INPUT.model_selection+1);
        set(handles.HTS_type,'Value',INPUT.HTS_type+1);

        %% Field Variables:
        
         Xfield(1) = INPUT.field_azimuth_or; % rectangular field shape azimuthial angle tilt from north [Deg]
         Xfield(2) = INPUT.field_elevation_or; % field plane elevation angle tilt from horizon @ azimuthial angle tilt from north [Deg]
         Xfield(3) = INPUT.N_dishes_per_cluster;%ceil(X0(3)); % number of dishes at each cluster (minimum 1)
         Xfield(4) = INPUT.N_clusters_in_field;%ceil(X0(4)); % number of clusters in the field (minimum 1)
         Xfield(5) = INPUT.alpha; % field parralelogram angle (from E-W line, +CW) [deg]
         Xfield(6) = INPUT.lns; % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
         Xfield(7) = INPUT.lew; % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
         Xfield(8) = INPUT.HEX_distance_from_field_center; % HEX distance from center field - HEX linkage [m]
         field_table(:,2) = num2cell(Xfield);
         set(handles.field_table,'Data',field_table);
        
        %% Dish & Receiver Variables:

        Xdish(1) = INPUT.Dish_effective_Area; % effective reflective area of each dish [m2]
        Xdish(2) = INPUT.Slope_Error;%sqrt(2.*2.5.^2); % mirror slope error (production quality) [mili.radians]
        Xdish(3) = INPUT.reflectivity.*100; % mirror max reflectivity [%]
        Xdish(4) = INPUT.receiver_peak_efficiency.*100; %[%]
         dishNreceiver_table(:,2) = num2cell(Xfield);
         set(handles.field_table,'Data',dishNreceiver_table);
        
        %% HTS Ducts Size (N1):
        XHTS = cell2mat(ducts_size_table(:,2));
        XHTS(1) = INPUT.HTS_hot.duct_R; % N1 hot pipe radius [m]
        XHTS(2) = INPUT.HTS_cold.duct_R; % N1 cold pipe radius [m]
         ducts_size_table(:,2) = num2cell(Xfield);
         set(handles.field_table,'Data',ducts_size_table);
 

        %% HTS Insulation Specification:
        
        if HTS_type>1
            Xthick(1:3) = INPUT.HTS_hot.insulation_occupation.*100;
            Xmaterial(1:3) = INPUT.HTS_hot.insulation_material;
            Xquality(3) = INPUT.HTS_hot.insulation_quality.*100;

            Xthick(4:6) = INPUT.HTS_cold.insulation_thickness;
            Xmaterial(4:6) = INPUT.HTS_cold.insulation_materia;
            Xquality(4:6) = INPUT.HTS_cold.insulation_quality.*100;
            
            annulus_insulation_table(:,2:4) = num2cell([Xthick,Xmaterial,Xquality]);
            insulation_table = annulus_insulation_table;
        else
            Xthick(1:3) = INPUT.HTS_hot.insulation_thickness;
            Xmaterial(1:3) = INPUT.HTS_hot.insulation_material;
            Xquality(3) = INPUT.HTS_hot.insulation_quality.*100;

            Xthick(4:6) = INPUT.HTS_cold.insulation_thickness;
            Xmaterial(4:6) = INPUT.HTS_cold.insulation_material;
            Xquality(4:6) = INPUT.HTS_cold.insulation_quality.*100;
            pipes_insulation_table(:,2:4) = num2cell([Xthick,Xmaterial,Xquality]);
            insulation_table = pipes_insulation_table;
        end
        set(handles.insulation_table,'Data',insulation_table);
        

                %% HEX Variables:
              
        XHEX(1,1) = INPUT.HP_pinch;
        XHEX(2,1) = INPUT.HP_approach;
        XHEX(3,1) = INPUT.HP_temperature;
        XHEX(4,1) = INPUT.HP_pressure;
        XHEX(1,2) = INPUT.MP_pinch;
        XHEX(2,2) = INPUT.MP_approach;
        XHEX(3,2) = INPUT.MP_temperature;
        XHEX(4,2) = INPUT.MP_pressure;
        XHEX(1,3) = INPUT.LP_pinch;
        XHEX(2,3) = INPUT.LP_approach;
        XHEX(3,3) = INPUT.LP_temperature;
        XHEX(4,3) = INPUT.LP_pressure;     
        HEX_table(:,2:4) = num2cell(XHEX);
        set(handles.HEX_table,'Data',HEX_table);

        
        %% Plant Variables:
        
        Xplant = cell2mat();
        Xplant(1) = INPUT.WS_pressure;
        Xplant(2) = INPUT.WS_temperature;
        Xplant(3) = INPUT.design_inlet_temperature;
        Xplant(4) = INPUT.power_block_peak.*100;
        plant_table(:,2) = num2cell(Xplant);
        set(handles.plant_table,'Data',plant_table);


        
if model_type>5
    

    %% Optimization Constraints:
    
     c(1) = INPUT.OPT_USER.A_max;
     b(1) = INPUT.OPT_USER.A_min;
     c(2) = INPUT.OPT_USER.B_max;
     b(2) = INPUT.OPT_USER.B_min;
     c(3) = INPUT.OPT_USER.C_max;
     b(3) = INPUT.OPT_USER.C_min;
     c(4) = INPUT.OPT_USER.D_max;
     b(4) = INPUT.OPT_USER.D_min;         
     c(5) = INPUT.OPT_USER.E_max;
     b(5) = INPUT.OPT_USER.E_min;
     c(6) = INPUT.OPT_USER.F_max;
     b(6) = INPUT.OPT_USER.F_min;
     c(7) = INPUT.OPT_USER.G_max;
     b(7) = INPUT.OPT_USER.G_min;
     c(8) = INPUT.OPT_USER.H_max;
     b(8) = INPUT.OPT_USER.H_min;  
    constraints_table(:,3:4) = [b',c'];
    set(handles.plant_table,'Data',constraints_table);
    
            
            %% Field Bounds Specification:

    num_field_Xb_table = [a',b',c'];
    
            c(1) = INPUT.UB.field_azimuth_or; % rectangular field shape azimuthial angle tilt from north [Deg]
            c(2) = INPUT.UB.field_elevation_or; % field plane elevation angle tilt from horizon @ azimuthial angle tilt from north [Deg]
            c(3) = INPUT.UB.N_dishes_per_cluster; % number of dishes at each cluster (minimum 1)
            c(4) = INPUT.UB.N_clusters_in_field; % number of clusters in the field (minimum 1)
            c(5) = INPUT.UB.alpha; % field parralelogram angle (from E-W line, +CW) [deg]
            c(6) = INPUT.UB.lns; % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
            c(7) = INPUT.UB.lew; % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
            
            b(1) = INPUT.LB.field_azimuth_or; % rectangular field shape azimuthial angle tilt from north [Deg]
            b(2) = INPUT.LB.field_elevation_or; % field plane elevation angle tilt from horizon @ azimuthial angle tilt from north [Deg]
            b(3) = INPUT.LB.N_dishes_per_cluster; % number of dishes at each cluster (minimum 1)
            b(4) = INPUT.LB.N_clusters_in_field; % number of clusters in the field (minimum 1)
            b(5) = INPUT.LB.alpha; % field parralelogram angle (from E-W line, +CW) [deg]
            b(6) = INPUT.LB.lns; % distance between dishes along field (if field is tilted, the axis tilts along) N-S axis [m]
            b(7) = INPUT.LB.lew; % distance between dishes along field (if field is tilted, the axis tilts along) E-W axis [m]
            field_Xb_table(:,3:4) = [b',c'];
            set(handles.HTS_Xb_table,'Data',field_Xb_table);
                        
            
            %% HTS Bounds Specification:

            c(1) = INPUT.UB.HTS_hot.duct_R; % N1 hot pipe radius [m]
            c(2) = INPUT.UB.HTS_cold.duct_R; % N1 cold pipe radius [m]
            c(3) = INPUT.UB.HTS_hot.insulation_occupation.*100;
            c(4) = INPUT.UB.HTS_hot.insulation_material;
            c(5) = INPUT.UB.HTS_hot.insulation_quality.*100;
            c(6) = INPUT.UB.HTS_cold.insulation_thickness;
            c(7) = INPUT.UB.HTS_cold.insulation_material;
            c(8) = INPUT.UB.HTS_cold.insulation_quality.*100;
            
            b(1) = INPUT.LB.HTS_hot.duct_R; % N1 hot pipe radius [m]
            b(1) = INPUT.LB.HTS_cold.duct_R; % N1 cold pipe radius [m]
            b(1) = INPUT.LB.HTS_hot.insulation_occupation.*100;
            b(1) = INPUT.LB.HTS_hot.insulation_material;
            b(1) = INPUT.LB.HTS_hot.insulation_quality;
            b(1) = INPUT.LB.HTS_cold.insulation_thickness;
            b(1) = INPUT.LB.HTS_cold.insulation_material;
            b(1) = INPUT.LB.HTS_cold.insulation_quality.*100;
            HTS_Xb_table(:,3:4) = [b',c'];
            set(handles.HTS_Xb_table,'Data',HTS_Xb_table);
            
            
end

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    HTS_type...
    file_name...
    error_count    

error_count = 0;
    check()
    
    if error_count==0; 
        harvest() 
        uisave
        refresh; refreshdata;
    end
    
       

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global field_table...
    dishNreceiver_table...
    ducts_size_table...
    annulus_insulation_table....
    pipes_insulation_table...
    insulation_table...
    HEX_table...
    plant_table... % optimization V
    field_Xb_table...
    annulus_Xb_table... 
    pipes_Xb_table...
    HTS_Xb_table...
    constraints_table...
    model_type...
    location...
    plant_type...
    system_pressure...
    HTS_type...
    file_name...
    error_count...
    stop_flag

    model_type = get(handles.model_type,'Value')-1;
    location = get(handles.location,'Value')-1;
    plant_type = get(handles.plant_type,'Value')-1;
    system_pressure = get(handles.system_pressure,'Value')-1;
    HTS_type = get(handles.HTS_type,'Value')-1;
    file_name = get(handles.file_name,'String');

    field_table = get(handles.field_table,'Data');
    dishNreceiver_table = get(handles.dishNreceiver_table,'Data');
    ducts_size_table = get(handles.ducts_size_table,'Data');
%     annulus_insulation_table = get(handles.annulus_insulation_table,'Data');
%     pipes_insulation_table = get(handles.pipes_insulation_table,'Data');
    insulation_table = get(handles.insulation_table,'Data');
    HEX_table = get(handles.HEX_table,'Data');
    plant_table = get(handles.plant_table,'Data');
        % optimization V
    field_Xb_table = get(handles.field_Xb_table,'Data');
%     annulus_Xb_table = get(handles.annulus_Xb_table,'Data');
%     pipes_Xb_table = get(handles.pipes_Xb_table,'Data');
    HTS_Xb_table = get(handles.HTS_Xb_table,'Data');
    constraints_table = get(handles.constraints_table,'Data');
    
%     if HTS_type>1
%         insulation_table = annulus_insulation_table;
%         HTS_Xb_table = annulus_Xb_table;
%     else
%         insulation_table = pipes_insulation_table;
%         HTS_Xb_table = pipes_Xb_table;
%     end
    
    stop_flag=0;
    error_count = 0;
    check()
    
    if error_count==0; 
        harvest()
        V4()
    end
    

% --- Executes on button press in pause.
function pause_Callback(hObject, eventdata, handles)
% hObject    handle to pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stat = get(hObject,'UserData');
if strcmp(stat,'Paused')
    set(hObject, 'UserData', 'Processing')
    set(hObject, 'String', 'Pause Model');
    uiresume;
%     message = ['Model Run Paused'];
%     set(handles.status,'String',message);
%     drawnow expose update
else
    set(hObject, 'UserData','Paused');
    set(hObject, 'String', 'Continue Model');
    uiwait;   
%     drawnow expose update
end
refresh; refreshdata;

% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function file_name_Callback(hObject, eventdata, handles)
% hObject    handle to file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of file_name as text
%        str2double(get(hObject,'String')) returns contents of file_name as a double
global stop_flag
stop_flag=1;
% handles = guidata(hObject);
% stop(handles.timer_handle);
% delete(handles.timer_handle);

stat = get(handles.pause,'UserData');
if strcmp(stat,'Paused')
    set(handles.pause,'UserData','Processing');
    set(handles.pause,'String','Stopped by User');
    set(handles.pause,'SelectionHighlight','off');
    uiresume;
end
refresh; refreshdata;

% --- Executes during object creation, after setting all properties.
function file_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%REPLACEINSELECTION Finds and replace strings in the selected text
%
%	replaceInSelection(oldSubstr, newSubstr)
%
% Replace all occurences of 'oldSubstr' with 'newSubstr' in the currently
% selected text in the MATLAB editor. Uses the MATLAB Editor API introduced
% in MATLAB 2011a and the built-in STRREP function.
%
%IN:
%	oldSubstr - the substring of text so search for 
%	newSubstr - the text to replace oldSubstr with
%
% See also STRREP

% Copyright (C) Sam Johnson 2012

function replaceInSelection(oldSubstr, newSubstr)

activeEditor = matlab.desktop.editor.getActive;
selectionPosition = activeEditor.Selection;

%convert the selection from Lines/Columns to index
startPos = matlab.desktop.editor.positionInLineToIndex(activeEditor, selectionPosition(1), selectionPosition(2));
endPos = matlab.desktop.editor.positionInLineToIndex(activeEditor, selectionPosition(3), selectionPosition(4));

newText = activeEditor.Text;

newTextPre = newText(1:startPos-1);
newTextPost = newText(endPos+1:end);

newTextRep = strrep(newText(startPos:min(endPos,numel(newText))), oldSubstr, newSubstr);

activeEditor.Text = [newTextPre newTextRep newTextPost];

%Re-select
[selectionPosition(3) selectionPosition(4)] = matlab.desktop.editor.indexToPositionInLine(activeEditor, endPos);
activeEditor.Selection = selectionPosition;

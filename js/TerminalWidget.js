WT_DECLARE_WT_MEMBER
(ButtonClickedJsSlot, Wt::JavaScriptFunction, "ButtonClickedJsSlot",
function( id, edit ){
  //Called when button is pressed, to allow you to proccess the text in JS, before
  //  uploading to server (e.g. c++ land)
  if( $(edit).hasClass("Wt-edit-emptyText") )
  {
    console.log( "No text entered." );
    return;
  }
  
  console.log( "will emit lineentered signal" );
  // Wt.emit( id, {name: 'lineentered'}, "JS proccessed string: " + edit.value.replace(/e/g, "eeee") + "\nanswer!" );
  Wt.emit( id, {name: 'lineentered'}, edit.value );
 $(id).scrollTop += $(id).scrollHeight;
}
);


WT_DECLARE_WT_MEMBER
(TerminalWidgetInit, Wt::JavaScriptFunction, "TerminalWidgetInit",
function( arg1 ){
    this.input_history_array = [];
    this.input_history_index = 0;
    this.saved_input = "";
    this.traversed_inputs = [];
}
);


WT_DECLARE_WT_MEMBER
(TerminalWidgetExampleClientJsFcn, Wt::JavaScriptFunction, "TerminalWidgetExampleClientJsFcn",
function( edit ){
 console.log("TerminalWidgetExampleClientJsFcn being called");
  if( $(edit).hasClass("Wt-edit-emptyText") )
  {
    console.log( "No text entered." );
    return;
  }
}
);

WT_DECLARE_WT_MEMBER
(TerminalWidgetSaveInput, Wt::JavaScriptFunction, "TerminalWidgetSaveInput",
     function( text ){
        if ( this.input_history_index >= this.input_history_array.length || this.input_history_array.length === 0 )
            this.saved_input = text;
        else if ( this.input_history_array[ this.input_history_index ] !== text )
            this.input_history_array[ this.input_history_index ] = text;
     }
 );

WT_DECLARE_WT_MEMBER
(TerminalWidgetAccessPreviousInput, Wt::JavaScriptFunction, "TerminalWidgetAccessPreviousInput",
    function( ){
        if ( this.input_history_array.length === 0 )
            return "{empty}";
 
         var index_already_saved = this.traversed_inputs.length > 0 && this.traversed_inputs[ this.traversed_inputs.length-1 ] < this.input_history_index;
         if ( !index_already_saved ) { this.traversed_inputs.push( [ this.input_history_index, this.input_history_array[ this.input_history_index ] ] ); }

        if ( this.input_history_index > 0 )
            return this.input_history_array[ --this.input_history_index ];
 
        return this.input_history_array[0];
    }
 );

WT_DECLARE_WT_MEMBER
(TerminalWidgetAccessNextInput, Wt::JavaScriptFunction, "TerminalWidgetAccessNextInput",
    function( ){
        if ( this.input_history_array.length === 0 && this.saved_input.length === 0 )
            return "";

        if ( !( this.input_history_index >= this.input_history_array.length ) ) {
            ++this.input_history_index;

            if ( this.input_history_index < this.input_history_array.length )
                return this.input_history_array[ this.input_history_index ];
        }
 
        return this.saved_input;
    }
 );

WT_DECLARE_WT_MEMBER
(TerminalWidgetChartClicked, Wt::JavaScriptFunction, "TerminalWidgetChartClicked",
function( edit, energy, count, pageX, pageY ){
    var text, cursor_start, cursor_end;                   // line edit variables
    var first_part, second_part, selected_text;            // text-manipulation variables
 
    /* Not quite sure what to do with pageX and pageY arguments. 
     For future reference we can use these arguments for a certain function.
     
     The default replacement for a selected text is the energy value clicked inside the chart.
     
     If the user is highlighting an argument for a bin value, then the bin value is replaced with the text.
     */
 
    edit.focus();
 
    text = edit.value;
    cursor_start = edit.selectionStart;
    cursor_end = edit.selectionEnd;
 
    first_part = text.slice( 0, edit.selectionStart );
    second_part = text.slice( edit.selectionEnd );
    selected_text = text.slice( edit.selectionStart, edit.selectionEnd );
 
    if ( selected_text.match( /energy/gi ) != null ) {
        edit.value = first_part +  energy.toString() + second_part;
        edit.selectionStart = first_part.length;
        edit.selectionEnd = first_part.length + energy.toString().length;
 
    } else if ( selected_text.match( /bin/gi ) != null ) {
        edit.value = first_part +  count.toString() + second_part;
        edit.selectionStart = first_part.length;
        edit.selectionEnd = first_part.length + count.toString().length;

    } else {
        edit.value = first_part + energy.toString() + second_part;
        edit.selectionStart = first_part.length;
        edit.selectionEnd = first_part.length +  energy.toString().length;
    }
}
);

WT_DECLARE_WT_MEMBER
(TerminalWidgetHandleEnterKey, Wt::JavaScriptFunction, "TerminalWidgetHandleEnterKey",
     function( text ){
         if ( text !== "{empty}" && text.trim().length > 0  ) {
 
             for ( var i = 0; i < this.traversed_inputs.length; i++ ) {
                if ( this.traversed_inputs[i][0] === this.input_history_index ) {
                     this.input_history_array[this.input_history_index] = this.traversed_inputs[i][1];
                     break;
                }
             }
 
             this.input_history_array.push( text );
             this.input_history_index = this.input_history_array.length;
 
             if ( this.saved_input.length > 0 )
                this.saved_input = "";
 
             this.traversed_inputs = [];
         }
     }
 );


WT_DECLARE_WT_MEMBER
(TerminalWidgetKeyPressed, Wt::JavaScriptFunction, "TerminalWidgetKeyPressed",
function( edit , event ){
 
    var text, cursor_start, cursor_end;     // line edit variables
    var first_part, second_part, selected_text;            // text-manipulation variables
 
    function areParenthesisBalanced( input ){
        if ( !input.includes('(') ) { return false; }
     
        var stack = [];
        function arePair( opening, closing ) { return opening === '(' && closing === ')';  }

        for ( var i = 0; i < input.length; i++ ) {
            if ( input.charAt( i ) === '(' ) {
                stack.push( input.charAt( i ) );

            } else if ( input.charAt( i ) === ')' ) {
                if ( stack.length === 0 || !arePair( stack[ stack.length-1 ], input.charAt( i ) ) ) {
                    return false;
                 } else { stack.pop(); }
            }
        }
         return stack.length === 0;
     }
 
    text = edit.value;
 
    cursor_start = edit.selectionStart;
    cursor_end = edit.selectionEnd;
 
    first_part = text.slice( 0, edit.selectionStart );
    second_part = text.slice( edit.selectionEnd );
    selected_text = text.slice( cursor_start, cursor_end );
 
     switch ( event.key ) {
         case '(':
             //   If the end of the cursor is either at the end or not in behind a letter, add auto-parentheses, else just add left parenthesis.
             event.preventDefault();
             if ( cursor_start < cursor_end ) {
                edit.value = first_part + '(' + selected_text + ')' + second_part;
                edit.selectionStart = first_part.length + 1;
                edit.selectionEnd = first_part.length + selected_text.length + 1;
                break;
             }
 
             var to_add = cursor_end >= text.length || !text.charAt( cursor_end ).match( /[a-z]/i ) || !text.charAt( cursor_end ) === ')' ? '()' : '(';
             edit.value = first_part + to_add + second_part;
             edit.selectionEnd = edit.selectionStart = cursor_start + 1;
             break;
         case '"':
             //   If the end of the cursor is either at the end or not in behind a letter, add auto-quotes, else just add a regular quote.
             event.preventDefault();
             if ( cursor_start < cursor_end ) {
                 edit.value = first_part + '"' + selected_text + '"' + second_part;
                 edit.selectionStart = first_part.length + 1;
                 edit.selectionEnd = first_part.length + selected_text.length + 1;
                 break;
             }
 
             var to_add = cursor_end >= text.length || !text.charAt( cursor_end ).match( /[a-z]/i ) ? '""' : '"  ';
             edit.value = first_part + to_add + second_part;
             edit.selectionEnd = edit.selectionStart = cursor_start + 1;
             break;
         case ')':                      // If parenthesis are empty and user is attempting to close parenthesis, just move cursor one step instead of adding another set.
             event.preventDefault();
             if ( /*!areParenthesisBalanced( text.substring( cursor_start-1, cursor_end+1 ) )*/
                 text.substring( cursor_start-1, cursor_end+1 ) !== '()') {      // right now only avoids inserting ')' when in '()'
                 edit.value = first_part + ')' + second_part;
             }
             edit.selectionEnd = edit.selectionStart = cursor_start + 1;
             break;
         case 'Enter': this.TerminalWidgetHandleEnterKey(text); break;
         default :     this.TerminalWidgetSaveInput( text );    break;
     }
}
);



WT_DECLARE_WT_MEMBER
(TerminalWidgetKeyDown, Wt::JavaScriptFunction, "TerminalWidgetKeyDown",
 function( edit , event ){
         var text = edit.value;
 
         var cursor_start = edit.selectionStart;
         var cursor_end = edit.selectionEnd;
 
         var first_part = text.slice( 0, edit.selectionStart );
         var second_part = text.slice( edit.selectionEnd+1 );
 
         
         switch ( event.key ) {
            case 'Backspace':
                if ( text.length > 0 && cursor_start === cursor_end &&                                       // Auto-delete closed phrases '()' and ' "" '
                    ( (text.charAt( cursor_start-1 ) === '(' && text.charAt( cursor_end ) === ')') ||
                      (text.charAt( cursor_start-1 ) === '"' && text.charAt( cursor_end ) === '"') )  ) {
                    edit.value = first_part + second_part;
                    edit.selectionEnd = edit.selectionStart = cursor_start;
                }
                break;
            case 'ArrowUp':
                this.TerminalWidgetSaveInput( text );
 
                var previous_input = this.TerminalWidgetAccessPreviousInput();
                if (previous_input !== "{empty}")
                    edit.value = previous_input;
 
                event.preventDefault();
                break;
            case 'ArrowDown' :
                this.TerminalWidgetSaveInput( text );
 
                var next_input = this.TerminalWidgetAccessNextInput();
                edit.value = next_input;
 
                event.preventDefault();
                break;
            default: break;
         }
 
 }
 );

WT_DECLARE_WT_MEMBER
(TerminalWidgetSelectFirstArgument, Wt::JavaScriptFunction, "TerminalWidgetSelectFirstArgument",
 function( command ){
    var index_of_left_parenthesis = command.indexOf('(');
    var index_of_right_parenthesis = command.indexOf(')');
 
    if ( index_of_left_parenthesis === -1 || index_of_right_parenthesis === -1 )
        return [ 0, command.length-1 ];
 
    var arguments = command.substr( index_of_left_parenthesis + 1 );
    if ( arguments.trim() === ')' )
        return [ command.length, command.length ];
 
 
    var start_index_first_argument = arguments.search(/[a-z]/i) + 1;
    var end_index_of_first_argument = arguments.trim().search(/\\s|,/) + start_index_first_argument - 1;
    return [ index_of_left_parenthesis + start_index_first_argument, index_of_left_parenthesis + end_index_of_first_argument + 1 ];
 }
 );


WT_DECLARE_WT_MEMBER
(TerminalWidgetCommandMenuItemSelected, Wt::JavaScriptFunction, "TerminalWidgetCommandMenuItemSelected",
 function( edit, command ){
     edit.focus();
     
     if( $(edit).hasClass("Wt-edit-emptyText") ) {
         edit.value = command;
         return;
     }
     
     var cursor_start = edit.selectionStart;
     var cursor_end = edit.selectionEnd;
     var text = edit.value;
     
     var first_part = text.slice( 0, edit.selectionStart );
     var second_part = text.slice( edit.selectionEnd );
     
     var new_cursor_positions = this.TerminalWidgetSelectFirstArgument( command );
     edit.value = first_part + command + second_part;
     edit.selectionStart = first_part.length + new_cursor_positions[0];
     edit.selectionEnd = first_part.length + new_cursor_positions[1];
 }
 );


WT_DECLARE_WT_MEMBER
(TerminalWidgetDarken, Wt::JavaScriptFunction, "TerminalWidgetDarken",
 function( enteredtxt, edit ){
        enteredtxt.style.backgroundColor = 'black';
        enteredtxt.style.color = 'white';
        enteredtxt.style.borderColor = 'white';
     
        edit.style.backgroundColor = 'black';
        edit.style.color = 'white';
 }
 );


WT_DECLARE_WT_MEMBER
(TerminalWidgetLighten, Wt::JavaScriptFunction, "TerminalWidgetLighten",
 function( enteredtxt, edit ){
         enteredtxt.style.backgroundColor = 'white';
         enteredtxt.style.color = 'black';
         enteredtxt.style.borderColor = 'black';
         
         edit.style.backgroundColor = 'white';
         edit.style.color = 'black';
 }
 );

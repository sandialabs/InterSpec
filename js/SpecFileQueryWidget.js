
WT_DECLARE_WT_MEMBER
(FileQueryInit, Wt::JavaScriptFunction, "FileQueryInit",
function( divid, parentid, additional_filters )
{
  var rules_basic = {
    condition: 'AND',
    rules: [{
      id: 'Manufacturer',
      operator: 'contains',
      size: 65,
      value: 'Acme'
    }]
  };
  
  var is_valid_str = function(value,rule){
    if( rule.operator.type === 'equal' || rule.operator.type === 'not equal')
      return true;
    
    if( !value || value.length < 1 )
      return "Search term can not be blank";
    
    if( rule.operator.type === 'regex' )
    {
      try{
        new RegExp(value);
        return true;
      }catch(e){
        return "Invalid regular expression: " + e.toString();
      }
    }
    
    return true;
  };
  
  
  var timeValObj = {
    format: /\\s*((\\s*\\+?\\s*((\\d+(\\.\\d*)?)|(\\.\\d*))\\s*(?:[Ee][+\\-]?\\d+)?\\s*\\s*(year|yr|y|day|d|hrs|hour|h|minutes|min|m|second|s|ms|microseconds|us|nanoseconds|ns)\\s*)|((\\+?\\d+:\\d\\d:\\d+(\\.\\d+)?)\\s*)\\s*)+/,
    messages: {
      format: 'Time durations can be specified like "3h 47m", "60s", "12:43:12.32", etc'
    }
  };
  
  var builder = $('#'+divid);
  
  builder.on('afterCreateRuleFilters.queryBuilder', function(rule,el) {
    //We get here after a rule is created (which includes when you add a group that creates a rule by default)
    //console.log( 'afterCreateRuleFilters.queryBuilder' );
    //console.log(rule);
    //console.log(el);
    
    el.$el.find('.rule-value-container input').on('change',function(){
      builder.queryBuilder('validate');
      console.log( 'validate on change' );
    } );
    
    el.$el.find('.rule-value-container input').on('keypress',function(){
      builder.queryBuilder('validate');
      console.log( 'keypress on change' );
    } );
  });
  
  
  var doUpdateVal = function(){
    var result = builder.queryBuilder('getRules',{ allow_invalid: true, skip_empty: true });
    var resultjson = (!$.isEmptyObject(result)) ? JSON.stringify(result, null, 2) : 'empty';
    var oldquery = builder.data('QueryJson');
    if( oldquery !== resultjson )
    {
      builder.data('QueryJson',resultjson);
      Wt.emit(parentid, 'fileSearchQueryChanged', resultjson);
    }
  };
  
  //
  builder.on('afterUpdateGroupCondition.queryBuilder afterUpdateRuleFilter.queryBuilder afterUpdateRuleOperator.queryBuilder afterUpdateRuleValue.queryBuilder', function(){
    //We get here after values are updated/changed, basically validate 'on the fly'
    doUpdateVal();
  });
  
  //Validate after adding/removing/moving groups
  builder.on('afterAddGroup.queryBuilder', function(event, group){
    //if( more than one (group+rule) )
    //   $('#'+group.id).find('.group-conditions').css('visibility','hidden');
    //else $('#'+group.id).find('.group-conditions').css('visibility','visible');
    //However, the left bar thing is still there...
    //console.log('afterAddGroup.queryBuilder');
    //console.log(group);
    
    doUpdateVal();
  });
  builder.on('afterDeleteGroup.queryBuilder', function(event, group){
    //if( more than one (group+rule) )
    //   $('#'+group.id).find('.group-conditions').css('visibility','hidden');
    //else $('#'+group.id).find('.group-conditions').css('visibility','visible');
    //console.log('afterDeleteGroup.queryBuilder');
    //console.log(group);
    doUpdateVal();
  });
  builder.on('afterMove.queryBuilder', function(event, node){ doUpdateVal(); });
  
  builder.on('afterDeleteRule.queryBuilder', function(event, group){
    doUpdateVal();
  });
  
  builder.on('getGroupTemplate.queryBuilder.filter', function(e, level) {
    //console.log( 'On filter: ' + e.value );
    //We get here when a group is added
  });
  
  builder.on('afterUpdateRuleFilter.queryBuilder', function(e, rule) {
    //We get here when field type is changed (e.g., from filename, to detector type).
    //console.log( 'afterUpdateRuleFilter.queryBuilder' );
    doUpdateVal();
  });
  
  builder.on('afterUpdateRuleOperator.queryBuilder', function(e, rule) {
    //We get here when the operator is updated (e.x., change from 'equal' to 'not equal')
    console.log( 'afterUpdateRuleOperator.queryBuilder' );
    doUpdateVal();
  });
  
  builder.on('afterCreateRuleOperators.queryBuilder', function(e) {
    //We get here when field type is changed (e.g., from filename, to detector type).
    console.log( 'afterCreateRuleOperators.queryBuilder' );
    doUpdateVal();
  });
  
  builder.on('afterCreateRuleInput.queryBuilder', function(e, rule) { /* Looks to be called when the widget is created. */ });
  builder.on('afterApplyRuleFlags.queryBuilder', function(e, rule) { /* We get here after nearly all operations. */ });
  builder.on('afterSetRules.queryBuilder', function(e, rule) { /* Looks to be called when the widget is created. */ });


  var select_filters = [{
    id: 'Filename',
    label: 'Filename',
    size: 65,
    type: 'string',
    placeholder: '(filename includes entire path to file)',
    operators: ['contains', 'equal', 'not equal', 'does not contain', 'begins with', 'does not begin with', 'ends with', 'does not end with', 'regex' ],
    validation: { callback: is_valid_str }
  },{
    id: 'Detector name',
    label: 'Detector name',
    type: 'string',
    placeholder: '(Detector name within detection system. Ex. "Aa1", "Ba2", etc.)',
    operators: ['contains', 'equal', 'not equal', 'does not contain', 'begins with', 'does not begin with', 'ends with', 'does not end with', 'regex' ],
    validation: { callback: is_valid_str }
  },{
    id: 'Serial number',
    label: 'Serial number',
    placeholder: '(serial number specified in file)',
    size: 65,
    type: 'string',
    operators: ['contains', 'equal', 'not equal', 'does not contain', 'begins with', 'does not begin with', 'ends with', 'does not end with', 'regex' ],
    validation: { callback: is_valid_str }
  },{
    id: 'Manufacturer',
    label: 'Manufacturer',
    size: 65,
    type: 'string',
    placeholder: '(Detector manufacturer, explicit or inferred; see also &#39;Detector system&#39;)',
    operators: ['contains', 'equal', 'not equal', 'does not contain', 'begins with', 'does not begin with', 'ends with', 'does not end with', 'regex' ],
    validation: { callback: is_valid_str }
  },{
    id: 'Model',
    label: 'Model',
    size: 65,
    placeholder: '(Detector model name, explicit or inferred; see also &#39;Detector system&#39;)',
    type: 'string',
    operators: ['contains', 'equal', 'not equal', 'does not contain', 'begins with', 'does not begin with', 'ends with', 'does not end with', 'regex' ],
    validation: { callback: is_valid_str }
  },{
    id: 'Detector system',
    label: 'Detector system',
    type: 'integer',
    input: 'select',
    values: {
      0: 'GR135',
      1: 'IdentiFINDER',
      2: 'IdentiFINDER-NG',
      3: 'IdentiFINDER-LaBr3',
      4: 'Detective - Unknown Model',
      5: 'Detective-EX',
      6: 'Detective-EX100',
      7: 'Detective-200',
      8: 'SAIC8',
      9: 'Falcon 5000',
      10: 'MicroDetective',
      11: 'MicroRaider',
      12: 'SAM940',
      13: 'SAM940LaBr3',
      14: 'SAM945',
      15: 'RS-701',
      16: 'RS-705',
      17: 'RSI - Other/Unspecified',
      18: 'RadHunterNaI',
      19: 'RadHunterLaBr3',
      20: 'RadEagle NaI 3x1',
      21: 'RadEagle CeBr3 2x1',
      22: 'RadEagle CeBr3 3x0.8',
      23: 'RadEagle LaBr3 2x1',
      24: 'SRPM-210',
      25: 'Detective X',
      26: 'Unknown'
    },
    operators: ['equal']
  },{
    id: 'Energy Cal Type',
    label: 'Energy Cal Type',
    type: 'integer',
    input: 'select',
    values: {
      0: 'Polynomial',
      1: 'Full Range Fraction',
      2: 'Lower Channel Edge',
      3: 'Unknown'
    },
    operators: ['equal']
  },{
    id: 'UUID',
    label: 'UUID',
    type: 'string',
    operators: ['contains', 'equal', 'not equal', 'does not contain', 'begins with', 'does not begin with', 'ends with', 'does not end with', 'regex' ],
    validation: { callback: is_valid_str }
  },{
    id: 'Remark',
    label: 'Remark',
    type: 'string',
    operators: ['contains', 'equal', 'not equal', 'does not contain', 'begins with', 'does not begin with', 'ends with', 'does not end with', 'regex' ],
    validation: { callback: is_valid_str }
  },{
    id: 'Location name',
    label: 'Location name',
    type: 'string',
    size: 65,
    placeholder: '(ex. "Lane 3", "Test Site", etc.)',
    operators: ['contains', 'equal', 'not equal', 'does not contain', 'begins with', 'does not begin with', 'ends with', 'does not end with', 'regex' ],
    validation: { callback: is_valid_str }
  },{
    id: 'Has RIID Analysis',
    label: 'Has RIID Analysis',
    type: 'integer',
    input: 'radio',
    values: {
      1: 'Yes',
      0: 'No'
    },
    default_value: 1,
    operators: ['equal']
  },{
    id: 'RIID Ana result',
    label: 'RIID Ana result',
    size: 65,
    placeholder: '(text anywhere within RIID results)',
    type: 'string',
    operators: ['contains', 'equal', 'not equal', 'does not contain', 'begins with', 'does not begin with', 'ends with', 'does not end with', 'regex' ],
    validation: { callback: is_valid_str }
  },{
    id: 'RIID IDed nuclide',
    label: 'RIID IDed nuclide',
    type: 'string',
    size: 65,
    placeholder: '(RIID nuclide results. Inputs of Co60, CO-60, 60C0, etc are equiv.)',
    operators: ['contains', 'does not contain' ],
    validation: { callback: is_valid_str }
  },{
    id: 'Aquisition mode',
    label: 'Aquisition mode',
    type: 'integer',
    input: 'radio',
    values: {
      0: 'Dwell',
      1: 'Portal or Search'
    },
    default_value: 0,
    operators: ['equal']
  },{
    id: 'Has Neutron',
    label: 'Has Neutron',
    type: 'integer',
    input: 'radio',
    values: {
      1: 'Yes',
      0: 'No'
    },
    default_value: 1,
    operators: ['equal']
  },{
    id: 'Neutron CPS',
    label: 'Neutron CPS',
    type: 'double',
    operators: ['greater','less','equal','not_equal']
  },{
    id: 'Gamma CPS',
    label: 'Gamma CPS',
    type: 'double',
    operators: ['greater','less','equal','not_equal']
  },{
    id: 'Has Dev. Pairs',
    label: 'Has Dev. Pairs',
    type: 'integer',
    input: 'radio',
    values: {
      1: 'Yes',
      0: 'No'
    },
    default_value: 1,
    operators: ['equal']
  },{
    id: 'Has GPS Info',
    label: 'Has GPS Info',
    type: 'integer',
    input: 'radio',
    values: {
      1: 'Yes',
      0: 'No'
    },
    default_value: 1,
    operators: ['equal']
  },{
    id: 'Sum live time',
    label: 'Sum live time',
    type: 'string',
    size: 65,
    placeholder: '(ex. "300s", "1h 5 min 8 sec", "01:14:16", etc)',
    validation: timeValObj,
    operators: ['greater','less','equal','not_equal']
  },{
    id: 'Sum real time',
    label: 'Sum real time',
    type: 'string',
    size: 65,
    placeholder: '(ex. "300s", "1h 5 min 8 sec", "01:14:16", etc)',
    validation: timeValObj,
    operators: ['greater','less','equal','not_equal']
  },{
    id: 'Spectrum live time',
    label: 'Spectrum live time',
    type: 'string',
    size: 65,
    placeholder: '(ex. "300s", "1h 5 min 8 sec", "01:14:16", etc)',
    validation: timeValObj,
    operators: ['greater','less','equal','not_equal']
  },{
    id: 'Spectrum real time',
    label: 'Spectrum real time',
    type: 'string',
    size: 65,
    placeholder: '(ex. "300s", "1h 5 min 8 sec", "01:14:16", etc)',
    validation: timeValObj,
    operators: ['greater','less','equal','not_equal']
  },{
    id: 'Start Time',
    label: 'Start Time',
    type: 'datetime',
    validation: {format:'YYYY/MM/DD HH:mm'},
    placeholder: 'YYYY/MM/DD HH:mm:ss',
    size: 65,
    validation: timeValObj,
    operators: ['equal','not_equal','less','greater']
  },{
    id: 'Num. Time Samples',
    label: 'Num. Time Samples',
    type: 'integer',
    size: 65,
    placeholder: '(Num time samples system recorded)',
    operators: ['equal','not_equal','less','greater']
  },{
    id: 'Num. Records',
    label: 'Num. Records',
    type: 'integer',
    size: 65,
    placeholder: '(Num records in spectrum file)',
    operators: ['equal','not_equal','less','greater']
  },{
    id: 'Num Gamma Channels',
    label: 'Num Gamma Channels',
    size: 65,
    type: 'integer',
    operators: ['greater','less','equal','not_equal']
  },{
    id: 'Max Gamma Energy',
    label: 'Max Gamma Energy',
    size: 65,
    placeholder: '(Energy in keV)',
    type: 'double',
    operators: ['greater','less','equal','not_equal']
  },{
    id: 'Latitude',
    label: 'Latitude',
    type: 'double',
    operators: ['greater','less','equal','not_equal']
  },{
    id: 'Longitude',
    label: 'Longitude',
    type: 'double',
    operators: ['greater','less','equal','not_equal']
  }];


  if( additional_filters.length )
    select_filters = select_filters.concat(additional_filters);
  
  builder.queryBuilder({
    //plugins: ['bt-tooltip-errors'],
    
    /*
    icons: {
      add_group: 'fas fa-plus-square',
      add_rule: 'fas fa-plus-circle',
      remove_group: 'fas fa-minus-square',
      remove_rule: 'fas fa-minus-circle',
      error: 'fas fa-exclamation-triangle'
    },
    */
    operators: $.fn.queryBuilder.constructor.DEFAULTS.operators.concat([
      { type: 'regex',  nb_inputs: 1, multiple: false, apply_to: ['string'] },
      { type: 'begins with',  nb_inputs: 1, multiple: false, apply_to: ['string'] },
      { type: 'ends with',  nb_inputs: 1, multiple: false, apply_to: ['string'] },
      { type: 'not equal',  nb_inputs: 1, multiple: false, apply_to: ['string'] },
      { type: 'does not contain',  nb_inputs: 1, multiple: false, apply_to: ['string'] },
      { type: 'does not begin with',  nb_inputs: 1, multiple: false, apply_to: ['string'] },
      { type: 'does not end with',  nb_inputs: 1, multiple: false, apply_to: ['string'] }/*,
      { type: 'distance from',  nb_inputs: 2, multiple: false, apply_to: ['string'] }*/
    ]),
    
    /*
     Each entry in the 'filters' arra must have an 'id' field that exactly chantches
     The string provided by #to_string(FileDataField).
    */
    
    filters: select_filters,
    
    rules: rules_basic
  });
  
  
  /*
  $('#btn-get').on('click', function() {
    var result = builder.queryBuilder('getRules');
    
    if (!$.isEmptyObject(result)) {
      alert(JSON.stringify(result, null, 2));
    }
  });
   */
}
);






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

  var builderEl = document.getElementById( divid );
  var qb = null;

  var doUpdateVal = function(){
    if( !qb )
      return;
    var result = qb.getRules( { allow_invalid: true, skip_empty: true } );
    var resultjson = (!result || Object.keys(result).length === 0) ? 'empty' : JSON.stringify(result, null, 2);
    var oldquery = builderEl._QueryJson;
    if( oldquery !== resultjson )
    {
      builderEl._QueryJson = resultjson;
      Wt.emit(parentid, 'fileSearchQueryChanged', resultjson);
    }
  };

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
      4: 'IdentiFINDER-Tungsten',
      5: 'IdentiFINDER-R500 NaI',
      6: 'IdentiFINDER-R500 LaBr2',
      7: 'IdentiFINDER Ambiguous',
      8: 'Detective - Ambiguous',
      9: 'Detective-EX/DX',
      10: 'Detective-EX100',
      11: 'Detective-200',
      12: 'Detective X',
      13: 'MicroDetective',
      14: 'SAIC8',
      15: 'Falcon 5000',
      16: 'MicroRaider',
      17: 'SAM940',
      18: 'SAM940LaBr3',
      19: 'SAM945',
      20: 'RS-701',
      21: 'RS-705',
      22: 'RSI - Other/Unspecified',
      23: 'RadHunterNaI',
      24: 'RadHunterLaBr3',
      25: 'Interceptor',
      26: 'RIID-Eye-NaI',
      27: 'RIID-Eye-LaBr3',
      28: 'RadSeeker-NaI',
      29: 'RadSeeker-LaBr3',
      30: 'RadEagle NaI 3x1',
      31: 'RadEagle CeBr3 2x1',
      32: 'RadEagle CeBr3 3x0.8',
      33: 'RadEagle LaBr3 2x1',
      34: 'SRPM-210',
      35: 'Verifinder-NaI',
      36: 'Verifinder-LaBr3',
      37: 'Kromek D3S',
      38: 'RadiaCode CsI 10X',
      39: 'RadiaCode CsI 110',
      40: 'RadiaCode 103G',
      41: 'Fulcrum',
      42: 'Fulcrum40h',
      43: 'IdentiFinder-R425-NaI',
      44: 'IdentiFinder-R425-LaBr',
      45: 'Sam-950',
      46: 'KromekD5',
      47: 'KromekGR1',
      48: 'Raysid',
      49: 'Unknown'
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
    type: 'string',
    placeholder: 'YYYY/MM/DD HH:mm:ss',
    size: 65,
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

  qb = QueryBuilder.create( builderEl, {
    operators: QueryBuilder.DEFAULTS.operators.concat([
      { type: 'regex',  nb_inputs: 1, multiple: false, apply_to: ['string'] },
      { type: 'begins with',  nb_inputs: 1, multiple: false, apply_to: ['string'] },
      { type: 'ends with',  nb_inputs: 1, multiple: false, apply_to: ['string'] },
      { type: 'not equal',  nb_inputs: 1, multiple: false, apply_to: ['string'] },
      { type: 'does not contain',  nb_inputs: 1, multiple: false, apply_to: ['string'] },
      { type: 'does not begin with',  nb_inputs: 1, multiple: false, apply_to: ['string'] },
      { type: 'does not end with',  nb_inputs: 1, multiple: false, apply_to: ['string'] }
    ]),

    filters: select_filters,
    rules: rules_basic,

    events: {
      'afterCreateRuleFilters': function(e, rule) {
        rule.el.querySelectorAll('.rule-value-container input').forEach(function(input){
          input.addEventListener('change', function(){ if(qb) qb.validate(); });
          input.addEventListener('keypress', function(){ if(qb) qb.validate(); });
        });
      },
      'afterUpdateGroupCondition': function(){ doUpdateVal(); },
      'afterUpdateRuleFilter': function(){ doUpdateVal(); },
      'afterUpdateRuleOperator': function(){ doUpdateVal(); },
      'afterUpdateRuleValue': function(){ doUpdateVal(); },
      'afterAddGroup': function(){ doUpdateVal(); },
      'afterDeleteGroup': function(){ doUpdateVal(); },
      'afterMove': function(){ doUpdateVal(); },
      'afterDeleteRule': function(){ doUpdateVal(); },
      'afterCreateRuleOperators': function(){ doUpdateVal(); }
    }
  });
}
);

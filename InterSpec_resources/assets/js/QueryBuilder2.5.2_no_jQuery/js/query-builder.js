/*!
 * QueryBuilder 2.5.2 (no jQuery)
 * Based on jQuery QueryBuilder 2.5.2 by Damien "Mistic" Sorel (http://www.strangeplanet.fr)
 * Licensed under MIT (https://opensource.org/licenses/MIT)
 *
 * Converted to vanilla JavaScript - no jQuery, doT.js, or other external dependencies.
 * Stripped to only the features used by InterSpec (no plugins, no optgroups, etc.).
 */
(function(root) {
"use strict";

// =========================================================================
// Utilities
// =========================================================================

/**
 * Checks if a value is a plain object (not an array, not null, prototype is Object)
 */
function isPlainObject( obj )
{
  return typeof obj === 'object' && obj !== null
      && Object.getPrototypeOf( obj ) === Object.prototype;
}

/**
 * Checks if an object has no own properties
 */
function isEmptyObject( obj )
{
  if( !obj || typeof obj !== 'object' )
    return true;
  return Object.keys( obj ).length === 0;
}

/**
 * Deep merge with 'replace' semantics for arrays.
 * Arrays in sources completely replace arrays in target rather than being merged.
 */
function deepMergeReplace( target )
{
  for( var i = 1; i < arguments.length; i++ )
  {
    var source = arguments[i];
    if( source === null || source === undefined )
      continue;

    if( Array.isArray( source ) )
    {
      // Replace mode: return a deep copy of the source array
      target = deepCopyArray( source );
      continue;
    }

    for( var key in source )
    {
      if( !source.hasOwnProperty( key ) )
        continue;

      var srcVal = source[key];
      if( target === srcVal )
        continue;

      if( srcVal !== undefined )
      {
        if( isPlainObject( srcVal ) )
        {
          var tgtVal = target[key];
          target[key] = deepMergeReplace( isPlainObject( tgtVal ) ? tgtVal : {}, srcVal );
        }
        else if( Array.isArray( srcVal ) )
        {
          target[key] = deepCopyArray( srcVal );
        }
        else
        {
          target[key] = srcVal;
        }
      }
    }
  }

  return target;
}

/**
 * Deep copies an array, recursively cloning objects and arrays within.
 */
function deepCopyArray( arr )
{
  var result = [];
  for( var i = 0; i < arr.length; i++ )
  {
    var item = arr[i];
    if( Array.isArray( item ) )
      result.push( deepCopyArray( item ) );
    else if( isPlainObject( item ) )
      result.push( deepMergeReplace( {}, item ) );
    else
      result.push( item );
  }
  return result;
}

/**
 * Shallow extend (like Object.assign but works in older browsers)
 */
function extend( target )
{
  for( var i = 1; i < arguments.length; i++ )
  {
    var source = arguments[i];
    if( source )
    {
      for( var key in source )
      {
        if( source.hasOwnProperty( key ) )
          target[key] = source[key];
      }
    }
  }
  return target;
}

/**
 * Deep extend (no replace semantics for arrays - arrays are deep copied)
 */
function deepExtend( target )
{
  for( var i = 1; i < arguments.length; i++ )
  {
    var source = arguments[i];
    if( !source )
      continue;

    for( var key in source )
    {
      if( !source.hasOwnProperty( key ) )
        continue;

      var srcVal = source[key];
      if( target === srcVal )
        continue;

      if( isPlainObject( srcVal ) )
      {
        var tgtVal = target[key];
        target[key] = deepExtend( isPlainObject( tgtVal ) ? tgtVal : {}, srcVal );
      }
      else if( Array.isArray( srcVal ) )
      {
        target[key] = deepCopyArray( srcVal );
      }
      else if( srcVal !== undefined )
      {
        target[key] = srcVal;
      }
    }
  }
  return target;
}

/**
 * Creates a DOM element from an HTML string
 */
function htmlToElement( html )
{
  var template = document.createElement( 'template' );
  template.innerHTML = html.trim();
  return template.content.firstChild;
}


// =========================================================================
// MiniEmitter - simple event emitter for the Model class
// =========================================================================

function MiniEmitter()
{
  this._handlers = {};
}

MiniEmitter.prototype.on = function( type, handler )
{
  // Support on({type: handler, ...}) syntax
  if( typeof type === 'object' )
  {
    for( var key in type )
    {
      if( type.hasOwnProperty( key ) )
        this.on( key, type[key] );
    }
    return this;
  }

  if( !this._handlers[type] )
    this._handlers[type] = [];
  this._handlers[type].push( handler );
  return this;
};

MiniEmitter.prototype.off = function( type, handler )
{
  if( !this._handlers[type] )
    return this;

  if( !handler )
  {
    delete this._handlers[type];
  }
  else
  {
    var handlers = this._handlers[type];
    for( var i = handlers.length - 1; i >= 0; i-- )
    {
      if( handlers[i] === handler )
        handlers.splice( i, 1 );
    }
  }
  return this;
};

MiniEmitter.prototype.trigger = function( type )
{
  var args = Array.prototype.slice.call( arguments, 1 );
  var handlers = this._handlers[type];
  if( handlers )
  {
    for( var i = 0; i < handlers.length; i++ )
      handlers[i].apply( null, args );
  }
};


// =========================================================================
// QueryBuilder Event Object
// =========================================================================

function QBEvent( type, builder )
{
  this.type = type;
  this.builder = builder;
  this._defaultPrevented = false;
}

QBEvent.prototype.preventDefault = function()
{
  this._defaultPrevented = true;
};

QBEvent.prototype.isDefaultPrevented = function()
{
  return this._defaultPrevented;
};


// =========================================================================
// Namespace
// =========================================================================
var Utils = {};


// =========================================================================
// QueryBuilder Constructor
// =========================================================================

/**
 * @param {HTMLElement} el
 * @param {object} options
 * @constructor
 */
var QueryBuilder = function( el, options )
{
  el.queryBuilder = this;

  /**
   * Container element (raw DOM)
   */
  this.el = el;

  /**
   * Configuration object
   */
  this.settings = deepMergeReplace( {}, QueryBuilder.DEFAULTS, options );

  /**
   * Internal model
   */
  this.model = new Model();

  /**
   * Internal status
   */
  this.status = {
    id: null,
    generated_id: false,
    group_id: 0,
    rule_id: 0
  };

  /**
   * List of filters
   */
  this.filters = this.settings.filters;

  /**
   * List of icons
   */
  this.icons = this.settings.icons;

  /**
   * List of operators
   */
  this.operators = this.settings.operators;

  /**
   * Templates
   */
  this.templates = this.settings.templates;

  /**
   * Translations object
   */
  this.lang = null;

  /**
   * Event handlers storage: { eventType: [handler, ...] }
   */
  this._events = {};

  /**
   * Stored references for DOM event listeners (for cleanup)
   */
  this._domListeners = [];

  // translations
  if( QueryBuilder.regional['en'] === undefined )
  {
    Utils.error( 'Config', '"en" translations not loaded.' );
  }
  this.lang = deepMergeReplace( {}, QueryBuilder.regional['en'],
    QueryBuilder.regional[this.settings.lang_code], this.settings.lang );

  // "allow_groups" can be boolean or int
  if( this.settings.allow_groups === false )
    this.settings.allow_groups = 0;
  else if( this.settings.allow_groups === true )
    this.settings.allow_groups = -1;

  // init templates - ensure all templates are functions
  var self = this;
  Object.keys( this.templates ).forEach( function( tpl ) {
    if( !self.templates[tpl] )
      self.templates[tpl] = QueryBuilder.templates[tpl];
  });

  // ensure we have a container id
  if( !this.el.id )
  {
    this.el.id = 'qb_' + Math.floor( Math.random() * 99999 );
    this.status.generated_id = true;
  }
  this.status.id = this.el.id;

  // INIT
  this.el.classList.add( 'query-builder' );
  this.el.classList.add( 'form-inline' );

  this.filters = this.checkFilters( this.filters );
  this.operators = this.checkOperators( this.operators );
  this.bindEvents();
};


// =========================================================================
// Event methods on QueryBuilder prototype
// =========================================================================

/**
 * Triggers an event on the builder
 */
QueryBuilder.prototype.trigger = function( type )
{
  var event = new QBEvent( this._toEventName( type ), this );
  var args = Array.prototype.slice.call( arguments, 1 );
  args.unshift( event );

  var handlers = this._events[event.type];
  if( handlers )
  {
    for( var i = 0; i < handlers.length; i++ )
      handlers[i].apply( null, args );
  }

  return event;
};

/**
 * Triggers an event and returns the modified value
 */
QueryBuilder.prototype.change = function( type, value )
{
  var event = new QBEvent( this._toEventName( type, true ), this );
  event.value = value;

  var args = Array.prototype.slice.call( arguments, 2 );
  args.unshift( event );

  var handlers = this._events[event.type];
  if( handlers )
  {
    for( var i = 0; i < handlers.length; i++ )
      handlers[i].apply( null, args );
  }

  return event.value;
};

/**
 * Attaches an event listener
 */
QueryBuilder.prototype.on = function( type, cb )
{
  var name = this._toEventName( type );
  if( !this._events[name] )
    this._events[name] = [];
  this._events[name].push( cb );
  return this;
};

/**
 * Removes an event listener
 */
QueryBuilder.prototype.off = function( type, cb )
{
  if( !type )
  {
    // Remove all
    this._events = {};
    return this;
  }

  var name = this._toEventName( type );
  if( !cb )
  {
    delete this._events[name];
  }
  else if( this._events[name] )
  {
    var handlers = this._events[name];
    for( var i = handlers.length - 1; i >= 0; i-- )
    {
      if( handlers[i] === cb )
        handlers.splice( i, 1 );
    }
  }
  return this;
};

/**
 * Appends '.queryBuilder' and optionally '.filter' to event names
 */
QueryBuilder.prototype._toEventName = function( name, filter )
{
  return name.split( ' ' ).map( function( type ) {
    return type + '.queryBuilder' + (filter ? '.filter' : '');
  }).join( ' ' );
};


// =========================================================================
// Static properties
// =========================================================================

QueryBuilder.types = {
  'string':   'string',
  'integer':  'number',
  'double':   'number',
  'date':     'datetime',
  'time':     'datetime',
  'datetime': 'datetime',
  'boolean':  'boolean'
};

QueryBuilder.inputs = [
  'text',
  'number',
  'radio',
  'checkbox',
  'select'
];

QueryBuilder.modifiable_options = [
  'display_errors',
  'allow_groups',
  'allow_empty',
  'default_condition',
  'default_filter'
];

QueryBuilder.selectors = {
  group_container:      '.rules-group-container',
  rule_container:       '.rule-container',
  filter_container:     '.rule-filter-container',
  operator_container:   '.rule-operator-container',
  value_container:      '.rule-value-container',
  error_container:      '.error-container',
  condition_container:  '.rules-group-header .group-conditions',

  rule_header:          '.rule-header',
  group_header:         '.rules-group-header',
  group_actions:        '.group-actions',
  rule_actions:         '.rule-actions',

  rules_list:           '.rules-group-body>.rules-list',

  group_condition:      '.rules-group-header [name$=_cond]',
  rule_filter:          '.rule-filter-container [name$=_filter]',
  rule_operator:        '.rule-operator-container [name$=_operator]',
  rule_value:           '.rule-value-container [name*=_value_]',

  add_rule:             '[data-add=rule]',
  delete_rule:          '[data-delete=rule]',
  add_group:            '[data-add=group]',
  delete_group:         '[data-delete=group]'
};

QueryBuilder.templates = {};
QueryBuilder.regional = {};

QueryBuilder.OPERATORS = {
  equal:            { type: 'equal',            nb_inputs: 1, multiple: false, apply_to: ['string', 'number', 'datetime', 'boolean'] },
  not_equal:        { type: 'not_equal',        nb_inputs: 1, multiple: false, apply_to: ['string', 'number', 'datetime', 'boolean'] },
  in:               { type: 'in',               nb_inputs: 1, multiple: true,  apply_to: ['string', 'number', 'datetime'] },
  not_in:           { type: 'not_in',           nb_inputs: 1, multiple: true,  apply_to: ['string', 'number', 'datetime'] },
  less:             { type: 'less',             nb_inputs: 1, multiple: false, apply_to: ['number', 'datetime'] },
  less_or_equal:    { type: 'less_or_equal',    nb_inputs: 1, multiple: false, apply_to: ['number', 'datetime'] },
  greater:          { type: 'greater',          nb_inputs: 1, multiple: false, apply_to: ['number', 'datetime'] },
  greater_or_equal: { type: 'greater_or_equal', nb_inputs: 1, multiple: false, apply_to: ['number', 'datetime'] },
  between:          { type: 'between',          nb_inputs: 2, multiple: false, apply_to: ['number', 'datetime'] },
  not_between:      { type: 'not_between',      nb_inputs: 2, multiple: false, apply_to: ['number', 'datetime'] },
  begins_with:      { type: 'begins_with',      nb_inputs: 1, multiple: false, apply_to: ['string'] },
  not_begins_with:  { type: 'not_begins_with',  nb_inputs: 1, multiple: false, apply_to: ['string'] },
  contains:         { type: 'contains',         nb_inputs: 1, multiple: false, apply_to: ['string'] },
  not_contains:     { type: 'not_contains',     nb_inputs: 1, multiple: false, apply_to: ['string'] },
  ends_with:        { type: 'ends_with',        nb_inputs: 1, multiple: false, apply_to: ['string'] },
  not_ends_with:    { type: 'not_ends_with',    nb_inputs: 1, multiple: false, apply_to: ['string'] },
  is_empty:         { type: 'is_empty',         nb_inputs: 0, multiple: false, apply_to: ['string'] },
  is_not_empty:     { type: 'is_not_empty',     nb_inputs: 0, multiple: false, apply_to: ['string'] },
  is_null:          { type: 'is_null',          nb_inputs: 0, multiple: false, apply_to: ['string', 'number', 'datetime', 'boolean'] },
  is_not_null:      { type: 'is_not_null',      nb_inputs: 0, multiple: false, apply_to: ['string', 'number', 'datetime', 'boolean'] }
};

QueryBuilder.DEFAULTS = {
  filters: [],

  sort_filters: false,
  display_errors: true,
  allow_groups: -1,
  allow_empty: false,
  conditions: ['AND', 'OR'],
  default_condition: 'AND',
  inputs_separator: ' , ',
  select_placeholder: '------',
  display_empty_filter: true,
  default_filter: null,

  default_rule_flags: {
    filter_readonly: false,
    operator_readonly: false,
    value_readonly: false,
    no_delete: false
  },

  default_group_flags: {
    condition_readonly: false,
    no_add_rule: false,
    no_add_group: false,
    no_delete: false
  },

  templates: {
    group: null,
    rule: null,
    filterSelect: null,
    operatorSelect: null,
    ruleValueSelect: null
  },

  lang_code: 'en',
  lang: {},

  operators: [
    'equal',
    'not_equal',
    'in',
    'not_in',
    'less',
    'less_or_equal',
    'greater',
    'greater_or_equal',
    'between',
    'not_between',
    'begins_with',
    'not_begins_with',
    'contains',
    'not_contains',
    'ends_with',
    'not_ends_with',
    'is_empty',
    'is_not_empty',
    'is_null',
    'is_not_null'
  ],

  icons: {
    add_group:    'glyphicon glyphicon-plus-sign',
    add_rule:     'glyphicon glyphicon-plus',
    remove_group: 'glyphicon glyphicon-remove',
    remove_rule:  'glyphicon glyphicon-remove',
    error:        'glyphicon glyphicon-warning-sign'
  }
};


/**
 * Gets or extends the default configuration
 */
QueryBuilder.defaults = function( options )
{
  if( typeof options == 'object' )
  {
    deepMergeReplace( QueryBuilder.DEFAULTS, options );
  }
  else if( typeof options == 'string' )
  {
    if( typeof QueryBuilder.DEFAULTS[options] == 'object' )
      return deepExtend( {}, QueryBuilder.DEFAULTS[options] );
    else
      return QueryBuilder.DEFAULTS[options];
  }
  else
  {
    return deepExtend( {}, QueryBuilder.DEFAULTS );
  }
};

/**
 * Adds new methods to QueryBuilder prototype
 */
QueryBuilder.extend = function( methods )
{
  extend( QueryBuilder.prototype, methods );
};


// =========================================================================
// QueryBuilder prototype methods
// =========================================================================

/**
 * Final initialisation of the builder
 */
QueryBuilder.prototype.init = function( rules )
{
  this.trigger( 'afterInit' );

  if( rules )
  {
    this.setRules( rules );
    delete this.settings.rules;
  }
  else
  {
    this.setRoot( true );
  }
};

/**
 * Checks the configuration of each filter
 */
QueryBuilder.prototype.checkFilters = function( filters )
{
  var definedFilters = [];

  if( !filters || filters.length === 0 )
    Utils.error( 'Config', 'Missing filters list' );

  filters.forEach( function( filter, i ) {
    if( !filter.id )
      Utils.error( 'Config', 'Missing filter {0} id', i );

    if( definedFilters.indexOf( filter.id ) != -1 )
      Utils.error( 'Config', 'Filter "{0}" already defined', filter.id );

    definedFilters.push( filter.id );

    if( !filter.type )
      filter.type = 'string';
    else if( !QueryBuilder.types[filter.type] )
      Utils.error( 'Config', 'Invalid type "{0}"', filter.type );

    if( !filter.input )
      filter.input = QueryBuilder.types[filter.type] === 'number' ? 'number' : 'text';
    else if( QueryBuilder.inputs.indexOf( filter.input ) == -1 )
      Utils.error( 'Config', 'Invalid input "{0}"', filter.input );

    if( filter.operators )
    {
      filter.operators.forEach( function( operator ) {
        if( typeof operator != 'string' )
          Utils.error( 'Config', 'Filter operators must be global operators types (string)' );
      });
    }

    if( !filter.field )
      filter.field = filter.id;
    if( !filter.label )
      filter.label = filter.field;

    // Handle select filter values
    if( filter.input === 'select' )
    {
      var cleanValues = [];
      Utils.iterateOptions( filter.values, function( value, label ) {
        cleanValues.push( { value: value, label: label } );
      });
      filter.values = cleanValues;

      if( filter.placeholder )
      {
        if( filter.placeholder_value === undefined )
          filter.placeholder_value = -1;
      }
    }
  }, this );

  if( this.settings.sort_filters )
  {
    var self = this;
    filters.sort( function( a, b ) {
      return self.translate( a.label ).localeCompare( self.translate( b.label ) );
    });
  }

  return filters;
};

/**
 * Checks the configuration of each operator
 */
QueryBuilder.prototype.checkOperators = function( operators )
{
  var definedOperators = [];

  operators.forEach( function( operator, i ) {
    if( typeof operator == 'string' )
    {
      if( !QueryBuilder.OPERATORS[operator] )
        Utils.error( 'Config', 'Unknown operator "{0}"', operator );

      operators[i] = operator = deepMergeReplace( {}, QueryBuilder.OPERATORS[operator] );
    }
    else
    {
      if( !operator.type )
        Utils.error( 'Config', 'Missing "type" for operator {0}', i );

      if( QueryBuilder.OPERATORS[operator.type] )
        operators[i] = operator = deepMergeReplace( {}, QueryBuilder.OPERATORS[operator.type], operator );

      if( operator.nb_inputs === undefined || operator.apply_to === undefined )
        Utils.error( 'Config', 'Missing "nb_inputs" and/or "apply_to" for operator "{0}"', operator.type );
    }

    if( definedOperators.indexOf( operator.type ) != -1 )
      Utils.error( 'Config', 'Operator "{0}" already defined', operator.type );

    definedOperators.push( operator.type );
  }, this );

  return operators;
};

/**
 * Adds all events listeners to the builder
 */
QueryBuilder.prototype.bindEvents = function()
{
  var self = this;
  var Selectors = QueryBuilder.selectors;

  // Delegated change handler - checks element name attributes to identify target
  var changeHandler = function( e )
  {
    var name = e.target.name || '';

    // group condition change (name ends with _cond)
    if( name.indexOf( '_cond' ) > -1 && e.target.checked )
    {
      var groupEl = e.target.closest( Selectors.group_container );
      if( groupEl )
        self.getModel( groupEl ).condition = e.target.value;
    }

    // rule filter change (name ends with _filter)
    if( name.indexOf( '_filter' ) > -1 )
    {
      var ruleFilterEl = e.target.closest( Selectors.rule_container );
      if( ruleFilterEl )
        self.getModel( ruleFilterEl ).filter = self.getFilterById( e.target.value );
    }

    // rule operator change (name ends with _operator)
    if( name.indexOf( '_operator' ) > -1 )
    {
      var ruleOpEl = e.target.closest( Selectors.rule_container );
      if( ruleOpEl )
        self.getModel( ruleOpEl ).operator = self.getOperatorByType( e.target.value );
    }
  };

  this.el.addEventListener( 'change', changeHandler );
  this._domListeners.push( { type: 'change', handler: changeHandler } );

  // Delegated click handler
  var clickHandler = function( e )
  {
    var target;

    // add rule button
    target = e.target.closest( Selectors.add_rule );
    if( target )
    {
      var groupEl = target.closest( Selectors.group_container );
      if( groupEl )
        self.addRule( self.getModel( groupEl ) );
      return;
    }

    // delete rule button
    target = e.target.closest( Selectors.delete_rule );
    if( target )
    {
      var ruleEl = target.closest( Selectors.rule_container );
      if( ruleEl )
        self.deleteRule( self.getModel( ruleEl ) );
      return;
    }

    if( self.settings.allow_groups !== 0 )
    {
      // add group button
      target = e.target.closest( Selectors.add_group );
      if( target )
      {
        var addGroupEl = target.closest( Selectors.group_container );
        if( addGroupEl )
          self.addGroup( self.getModel( addGroupEl ) );
        return;
      }

      // delete group button
      target = e.target.closest( Selectors.delete_group );
      if( target )
      {
        var delGroupEl = target.closest( Selectors.group_container );
        if( delGroupEl )
          self.deleteGroup( self.getModel( delGroupEl ) );
        return;
      }
    }
  };

  this.el.addEventListener( 'click', clickHandler );
  this._domListeners.push( { type: 'click', handler: clickHandler } );

  // model events
  this.model.on({
    'drop': function( node ) {
      node.el.remove();
      self.refreshGroupsConditions();
    },
    'add': function( parent, node, index ) {
      var rulesList = parent.el.querySelector( ':scope > ' + QueryBuilder.selectors.rules_list );
      if( index === 0 )
      {
        rulesList.prepend( node.el );
      }
      else
      {
        var prevEl = parent.rules[index - 1].el;
        prevEl.after( node.el );
      }
      self.refreshGroupsConditions();
    },
    'move': function( node, group, index ) {
      node.el.remove();
      var rulesList = group.el.querySelector( ':scope > ' + QueryBuilder.selectors.rules_list );
      if( index === 0 )
      {
        rulesList.prepend( node.el );
      }
      else
      {
        var prevEl = group.rules[index - 1].el;
        prevEl.after( node.el );
      }
      self.refreshGroupsConditions();
    },
    'update': function( node, field, value, oldValue ) {
      if( node instanceof Rule )
      {
        switch( field )
        {
          case 'error':    self.updateError( node ); break;
          case 'flags':    self.applyRuleFlags( node ); break;
          case 'filter':   self.updateRuleFilter( node, oldValue ); break;
          case 'operator': self.updateRuleOperator( node, oldValue ); break;
          case 'value':    self.updateRuleValue( node, oldValue ); break;
        }
      }
      else
      {
        switch( field )
        {
          case 'error':     self.updateError( node ); break;
          case 'flags':     self.applyGroupFlags( node ); break;
          case 'condition': self.updateGroupCondition( node, oldValue ); break;
        }
      }
    }
  });
};

/**
 * Creates the root group
 */
QueryBuilder.prototype.setRoot = function( addRule, data, flags )
{
  addRule = (addRule === undefined || addRule === true);

  var group_id = this.nextGroupId();
  var groupEl = htmlToElement( this.getGroupTemplate( group_id, 1 ) );

  this.el.appendChild( groupEl );
  this.model.root = new Group( null, groupEl );
  this.model.root.model = this.model;

  this.model.root.data = data;
  this.model.root.flags = extend( {}, this.settings.default_group_flags, flags );
  this.model.root.condition = this.settings.default_condition;

  this.trigger( 'afterAddGroup', this.model.root );

  if( addRule )
    this.addRule( this.model.root );

  return this.model.root;
};

/**
 * Adds a new group
 */
QueryBuilder.prototype.addGroup = function( parent, addRule, data, flags )
{
  addRule = (addRule === undefined || addRule === true);

  var level = parent.level + 1;

  var e = this.trigger( 'beforeAddGroup', parent, addRule, level );
  if( e.isDefaultPrevented() )
    return null;

  var group_id = this.nextGroupId();
  var groupEl = htmlToElement( this.getGroupTemplate( group_id, level ) );
  var model = parent.addGroup( groupEl );

  model.data = data;
  model.flags = extend( {}, this.settings.default_group_flags, flags );
  model.condition = this.settings.default_condition;

  this.trigger( 'afterAddGroup', model );
  this.trigger( 'rulesChanged' );

  if( addRule )
    this.addRule( model );

  return model;
};

/**
 * Tries to delete a group
 */
QueryBuilder.prototype.deleteGroup = function( group )
{
  if( group.isRoot() )
    return false;

  var e = this.trigger( 'beforeDeleteGroup', group );
  if( e.isDefaultPrevented() )
    return false;

  var del = true;
  group.each( 'reverse', function( rule ) {
    del &= this.deleteRule( rule );
  }, function( group ) {
    del &= this.deleteGroup( group );
  }, this );

  if( del )
  {
    group.drop();
    this.trigger( 'afterDeleteGroup' );
    this.trigger( 'rulesChanged' );
  }

  return del;
};

/**
 * Performs actions when a group's condition changes
 */
QueryBuilder.prototype.updateGroupCondition = function( group, previousCondition )
{
  var condInputs = group.el.querySelectorAll( ':scope > ' + QueryBuilder.selectors.group_header + ' [name$=_cond]' );
  condInputs.forEach( function( input ) {
    input.checked = (input.value === group.condition);
    input.parentNode.classList.toggle( 'active', input.value === group.condition );
  });

  this.trigger( 'afterUpdateGroupCondition', group, previousCondition );
  this.trigger( 'rulesChanged' );
};

/**
 * Updates the visibility of conditions based on number of rules inside each group
 */
QueryBuilder.prototype.refreshGroupsConditions = function()
{
  var Selectors = QueryBuilder.selectors;
  (function walk( group ) {
    if( !group.flags || (group.flags && !group.flags.condition_readonly) )
    {
      var condInputs = group.el.querySelectorAll( ':scope > ' + Selectors.group_header + ' [name$=_cond]' );
      condInputs.forEach( function( input ) {
        input.disabled = group.rules.length <= 1;
        input.parentNode.classList.toggle( 'disabled', group.rules.length <= 1 );
      });
    }

    group.each( null, function( group ) {
      walk( group );
    }, this );
  }( this.model.root ));
};

/**
 * Adds a new rule
 */
QueryBuilder.prototype.addRule = function( parent, data, flags )
{
  var e = this.trigger( 'beforeAddRule', parent );
  if( e.isDefaultPrevented() )
    return null;

  var rule_id = this.nextRuleId();
  var ruleEl = htmlToElement( this.getRuleTemplate( rule_id ) );
  var model = parent.addRule( ruleEl );

  model.data = data;
  model.flags = extend( {}, this.settings.default_rule_flags, flags );

  this.trigger( 'afterAddRule', model );
  this.trigger( 'rulesChanged' );

  this.createRuleFilters( model );

  if( this.settings.default_filter || !this.settings.display_empty_filter )
  {
    model.filter = this.change( 'getDefaultFilter',
      this.getFilterById( this.settings.default_filter || this.filters[0].id ),
      model
    );
  }

  return model;
};

/**
 * Tries to delete a rule
 */
QueryBuilder.prototype.deleteRule = function( rule )
{
  if( rule.flags.no_delete )
    return false;

  var e = this.trigger( 'beforeDeleteRule', rule );
  if( e.isDefaultPrevented() )
    return false;

  rule.drop();
  this.trigger( 'afterDeleteRule' );
  this.trigger( 'rulesChanged' );

  return true;
};

/**
 * Creates the filters for a rule
 */
QueryBuilder.prototype.createRuleFilters = function( rule )
{
  var filters = this.change( 'getRuleFilters', this.filters, rule );
  var filterSelectHTML = this.getRuleFilterSelect( rule, filters );

  rule.el.querySelector( QueryBuilder.selectors.filter_container ).innerHTML = filterSelectHTML;

  this.trigger( 'afterCreateRuleFilters', rule );
  this.applyRuleFlags( rule );
};

/**
 * Creates the operators for a rule and init the rule operator
 */
QueryBuilder.prototype.createRuleOperators = function( rule )
{
  var operatorContainer = rule.el.querySelector( QueryBuilder.selectors.operator_container );
  operatorContainer.innerHTML = '';

  if( !rule.filter )
    return;

  var operators = this.getOperators( rule.filter );
  var operatorSelectHTML = this.getRuleOperatorSelect( rule, operators );

  operatorContainer.innerHTML = operatorSelectHTML;

  // set the operator without triggering update event
  if( rule.filter.default_operator )
    rule.__.operator = this.getOperatorByType( rule.filter.default_operator );
  else
    rule.__.operator = operators[0];

  var opSelect = rule.el.querySelector( QueryBuilder.selectors.rule_operator );
  if( opSelect )
    opSelect.value = rule.operator.type;

  this.trigger( 'afterCreateRuleOperators', rule, operators );
  this.applyRuleFlags( rule );
};

/**
 * Creates the main input for a rule
 */
QueryBuilder.prototype.createRuleInput = function( rule )
{
  var valueContainer = rule.el.querySelector( QueryBuilder.selectors.value_container );
  valueContainer.innerHTML = '';

  rule.__.value = undefined;

  if( !rule.filter || !rule.operator || rule.operator.nb_inputs === 0 )
    return;

  var self = this;
  var inputs = [];
  var filter = rule.filter;

  for( var i = 0; i < rule.operator.nb_inputs; i++ )
  {
    var ruleInputHTML = this.getRuleInput( rule, i );
    if( i > 0 )
    {
      var sep = document.createTextNode( this.settings.inputs_separator );
      valueContainer.appendChild( sep );
    }

    var ruleInputEl = htmlToElement( ruleInputHTML );
    valueContainer.appendChild( ruleInputEl );
    inputs.push( ruleInputEl );
  }

  valueContainer.style.display = '';

  inputs.forEach( function( input ) {
    var eventTypes = 'change';
    if( filter.input_event )
      eventTypes += ' ' + filter.input_event;

    eventTypes.split( ' ' ).forEach( function( evt ) {
      if( evt.trim() )
      {
        input.addEventListener( evt.trim(), function() {
          if( !rule._updating_input )
          {
            rule._updating_value = true;
            rule.value = self.getRuleInputValue( rule );
            rule._updating_value = false;
          }
        });
      }
    });
  });

  this.trigger( 'afterCreateRuleInput', rule );

  if( filter.default_value !== undefined )
  {
    rule.value = filter.default_value;
  }
  else
  {
    rule._updating_value = true;
    rule.value = self.getRuleInputValue( rule );
    rule._updating_value = false;
  }

  this.applyRuleFlags( rule );
};

/**
 * Performs action when a rule's filter changes
 */
QueryBuilder.prototype.updateRuleFilter = function( rule, previousFilter )
{
  this.createRuleOperators( rule );
  this.createRuleInput( rule );

  var filterSelect = rule.el.querySelector( QueryBuilder.selectors.rule_filter );
  if( filterSelect )
    filterSelect.value = rule.filter ? rule.filter.id : '-1';

  if( previousFilter && rule.filter && previousFilter.id !== rule.filter.id )
    rule.data = undefined;

  this.trigger( 'afterUpdateRuleFilter', rule, previousFilter );
  this.trigger( 'rulesChanged' );
};

/**
 * Performs actions when a rule's operator changes
 */
QueryBuilder.prototype.updateRuleOperator = function( rule, previousOperator )
{
  var valueContainer = rule.el.querySelector( QueryBuilder.selectors.value_container );

  if( !rule.operator || rule.operator.nb_inputs === 0 )
  {
    valueContainer.style.display = 'none';
    rule.__.value = undefined;
  }
  else
  {
    valueContainer.style.display = '';

    if( valueContainer.children.length === 0 || !previousOperator ||
        rule.operator.nb_inputs !== previousOperator.nb_inputs )
    {
      this.createRuleInput( rule );
    }
  }

  if( rule.operator )
  {
    var opSelect = rule.el.querySelector( QueryBuilder.selectors.rule_operator );
    if( opSelect )
      opSelect.value = rule.operator.type;

    rule.__.value = this.getRuleInputValue( rule );
  }

  this.trigger( 'afterUpdateRuleOperator', rule, previousOperator );
  this.trigger( 'rulesChanged' );
};

/**
 * Performs actions when rule's value changes
 */
QueryBuilder.prototype.updateRuleValue = function( rule, previousValue )
{
  if( !rule._updating_value )
    this.setRuleInputValue( rule, rule.value );

  this.trigger( 'afterUpdateRuleValue', rule, previousValue );
  this.trigger( 'rulesChanged' );
};

/**
 * Changes a rule's properties depending on its flags
 */
QueryBuilder.prototype.applyRuleFlags = function( rule )
{
  var flags = rule.flags;
  var Selectors = QueryBuilder.selectors;

  var filterEl = rule.el.querySelector( Selectors.rule_filter );
  if( filterEl ) filterEl.disabled = flags.filter_readonly;

  var operatorEl = rule.el.querySelector( Selectors.rule_operator );
  if( operatorEl ) operatorEl.disabled = flags.operator_readonly;

  rule.el.querySelectorAll( Selectors.rule_value ).forEach( function( el ) {
    el.disabled = flags.value_readonly;
  });

  if( flags.no_delete )
  {
    var delBtn = rule.el.querySelector( Selectors.delete_rule );
    if( delBtn ) delBtn.remove();
  }

  this.trigger( 'afterApplyRuleFlags', rule );
};

/**
 * Changes group's properties depending on its flags
 */
QueryBuilder.prototype.applyGroupFlags = function( group )
{
  var flags = group.flags;
  var Selectors = QueryBuilder.selectors;

  group.el.querySelectorAll( ':scope > ' + Selectors.group_header + ' [name$=_cond]' ).forEach( function( input ) {
    input.disabled = flags.condition_readonly;
    input.parentNode.classList.toggle( 'readonly', flags.condition_readonly );
  });

  if( flags.no_add_rule )
  {
    var addRuleBtn = group.el.querySelector( Selectors.add_rule );
    if( addRuleBtn ) addRuleBtn.remove();
  }
  if( flags.no_add_group )
  {
    var addGroupBtn = group.el.querySelector( Selectors.add_group );
    if( addGroupBtn ) addGroupBtn.remove();
  }
  if( flags.no_delete )
  {
    var delGroupBtn = group.el.querySelector( Selectors.delete_group );
    if( delGroupBtn ) delGroupBtn.remove();
  }

  this.trigger( 'afterApplyGroupFlags', group );
};

/**
 * Clears all errors markers
 */
QueryBuilder.prototype.clearErrors = function( node )
{
  node = node || this.model.root;
  if( !node ) return;

  node.error = null;

  if( node instanceof Group )
  {
    node.each( function( rule ) {
      rule.error = null;
    }, function( group ) {
      this.clearErrors( group );
    }, this );
  }
};

/**
 * Adds/Removes error on a Rule or Group
 */
QueryBuilder.prototype.updateError = function( node )
{
  if( this.settings.display_errors )
  {
    if( node.error === null )
    {
      node.el.classList.remove( 'has-error' );
    }
    else
    {
      var errorMessage = this.translate( 'errors', node.error[0] );
      errorMessage = Utils.fmt( errorMessage, node.error.slice( 1 ) );

      errorMessage = this.change( 'displayError', errorMessage, node.error, node );

      node.el.classList.add( 'has-error' );
      var errContainer = node.el.querySelector( QueryBuilder.selectors.error_container );
      if( errContainer )
        errContainer.setAttribute( 'title', errorMessage );
    }
  }
};

/**
 * Triggers a validation error event
 */
QueryBuilder.prototype.triggerValidationError = function( node, error, value )
{
  if( !Array.isArray( error ) )
    error = [error];

  var e = this.trigger( 'validationError', node, error, value );
  if( !e.isDefaultPrevented() )
    node.error = error;
};

/**
 * Destroys the builder
 */
QueryBuilder.prototype.destroy = function()
{
  this.trigger( 'beforeDestroy' );

  if( this.status.generated_id )
    this.el.removeAttribute( 'id' );

  this.clear();
  this.model = null;

  // Remove DOM event listeners
  var self = this;
  this._domListeners.forEach( function( listener ) {
    self.el.removeEventListener( listener.type, listener.handler );
  });
  this._domListeners = [];

  this._events = {};
  this.el.classList.remove( 'query-builder' );

  delete this.el._queryBuilder;
  delete this.el.queryBuilder;
};

/**
 * Clear all rules and resets the root group
 */
QueryBuilder.prototype.reset = function()
{
  var e = this.trigger( 'beforeReset' );
  if( e.isDefaultPrevented() )
    return;

  this.status.group_id = 1;
  this.status.rule_id = 0;

  this.model.root.empty();

  this.model.root.data = undefined;
  this.model.root.flags = extend( {}, this.settings.default_group_flags );
  this.model.root.condition = this.settings.default_condition;

  this.addRule( this.model.root );

  this.trigger( 'afterReset' );
  this.trigger( 'rulesChanged' );
};

/**
 * Clears all rules and removes the root group
 */
QueryBuilder.prototype.clear = function()
{
  var e = this.trigger( 'beforeClear' );
  if( e.isDefaultPrevented() )
    return;

  this.status.group_id = 0;
  this.status.rule_id = 0;

  if( this.model.root )
  {
    this.model.root.drop();
    this.model.root = null;
  }

  this.trigger( 'afterClear' );
  this.trigger( 'rulesChanged' );
};

/**
 * Modifies the builder configuration
 */
QueryBuilder.prototype.setOptions = function( options )
{
  var self = this;
  Object.keys( options ).forEach( function( opt ) {
    if( QueryBuilder.modifiable_options.indexOf( opt ) !== -1 )
      self.settings[opt] = options[opt];
  });
};

/**
 * Returns the model associated to a DOM object, or the root model
 */
QueryBuilder.prototype.getModel = function( target )
{
  if( !target )
    return this.model.root;
  else if( target instanceof Node )
    return target;
  else
    return target._queryBuilderModel;
};

/**
 * Validates the whole builder
 */
QueryBuilder.prototype.validate = function( options )
{
  options = extend( { skip_empty: false }, options );

  this.clearErrors();
  var self = this;

  var valid = (function parse( group ) {
    var done = 0;
    var errors = 0;

    group.each( function( rule ) {
      if( !rule.filter && options.skip_empty )
        return;

      if( !rule.filter )
      {
        self.triggerValidationError( rule, 'no_filter', null );
        errors++;
        return;
      }

      if( !rule.operator )
      {
        self.triggerValidationError( rule, 'no_operator', null );
        errors++;
        return;
      }

      if( rule.operator.nb_inputs !== 0 )
      {
        var valid = self.validateValue( rule, rule.value );
        if( valid !== true )
        {
          self.triggerValidationError( rule, valid, rule.value );
          errors++;
          return;
        }
      }

      done++;

    }, function( group ) {
      var res = parse( group );
      if( res === true ) done++;
      else if( res === false ) errors++;
    });

    if( errors > 0 )
      return false;
    else if( done === 0 && !group.isRoot() && options.skip_empty )
      return null;
    else if( done === 0 && (!self.settings.allow_empty || !group.isRoot()) )
    {
      self.triggerValidationError( group, 'empty_group', null );
      return false;
    }

    return true;

  }( this.model.root ));

  return this.change( 'validate', valid );
};

/**
 * Gets an object representing current rules
 */
QueryBuilder.prototype.getRules = function( options )
{
  options = extend( {
    get_flags: false,
    allow_invalid: false,
    skip_empty: false
  }, options );

  var valid = this.validate( options );
  if( !valid && !options.allow_invalid )
    return null;

  var self = this;

  var out = (function parse( group ) {
    var groupData = {
      condition: group.condition,
      rules: []
    };

    if( group.data )
      groupData.data = deepMergeReplace( {}, group.data );

    if( options.get_flags )
    {
      var flags = self.getGroupFlags( group.flags, options.get_flags === 'all' );
      if( !isEmptyObject( flags ) )
        groupData.flags = flags;
    }

    group.each( function( rule ) {
      if( !rule.filter && options.skip_empty )
        return;

      var value = null;
      if( !rule.operator || rule.operator.nb_inputs !== 0 )
        value = rule.value;

      var ruleData = {
        id: rule.filter ? rule.filter.id : null,
        field: rule.filter ? rule.filter.field : null,
        type: rule.filter ? rule.filter.type : null,
        input: rule.filter ? rule.filter.input : null,
        operator: rule.operator ? rule.operator.type : null,
        value: value
      };

      if( (rule.filter && rule.filter.data) || rule.data )
        ruleData.data = deepMergeReplace( {}, rule.filter.data, rule.data );

      if( options.get_flags )
      {
        var flags = self.getRuleFlags( rule.flags, options.get_flags === 'all' );
        if( !isEmptyObject( flags ) )
          ruleData.flags = flags;
      }

      groupData.rules.push( self.change( 'ruleToJson', ruleData, rule ) );

    }, function( model ) {
      var data = parse( model );
      if( data.rules.length !== 0 || !options.skip_empty )
        groupData.rules.push( data );
    }, this );

    return self.change( 'groupToJson', groupData, group );

  }( this.model.root ));

  out.valid = valid;

  return this.change( 'getRules', out );
};

/**
 * Sets rules from object
 */
QueryBuilder.prototype.setRules = function( data, options )
{
  options = extend( { allow_invalid: false }, options );

  if( Array.isArray( data ) )
  {
    data = {
      condition: this.settings.default_condition,
      rules: data
    };
  }

  if( !data || !data.rules || (data.rules.length === 0 && !this.settings.allow_empty) )
    Utils.error( 'RulesParse', 'Incorrect data object passed' );

  this.clear();
  this.setRoot( false, data.data, this.parseGroupFlags( data ) );

  data = this.change( 'setRules', data, options );

  var self = this;

  (function add( data, group ) {
    if( group === null )
      return;

    if( data.condition === undefined )
      data.condition = self.settings.default_condition;
    else if( self.settings.conditions.indexOf( data.condition ) == -1 )
    {
      Utils.error( !options.allow_invalid, 'UndefinedCondition', 'Invalid condition "{0}"', data.condition );
      data.condition = self.settings.default_condition;
    }

    group.condition = data.condition;

    data.rules.forEach( function( item ) {
      var model;

      if( item.rules !== undefined )
      {
        if( self.settings.allow_groups !== -1 && self.settings.allow_groups < group.level )
        {
          Utils.error( !options.allow_invalid, 'RulesParse', 'No more than {0} groups are allowed', self.settings.allow_groups );
          self.reset();
        }
        else
        {
          model = self.addGroup( group, false, item.data, self.parseGroupFlags( item ) );
          if( model === null )
            return;
          add( item, model );
        }
      }
      else
      {
        if( !item.empty )
        {
          if( item.id === undefined )
          {
            Utils.error( !options.allow_invalid, 'RulesParse', 'Missing rule field id' );
            item.empty = true;
          }
          if( item.operator === undefined )
            item.operator = 'equal';
        }

        model = self.addRule( group, item.data, self.parseRuleFlags( item ) );
        if( model === null )
          return;

        if( !item.empty )
          model.filter = self.getFilterById( item.id, !options.allow_invalid );

        if( model.filter )
        {
          model.operator = self.getOperatorByType( item.operator, !options.allow_invalid );

          if( !model.operator )
            model.operator = self.getOperators( model.filter )[0];
        }

        if( model.operator && model.operator.nb_inputs !== 0 )
        {
          if( item.value !== undefined )
            model.value = item.value;
          else if( model.filter.default_value !== undefined )
            model.value = model.filter.default_value;
        }

        if( self.change( 'jsonToRule', model, item ) != model )
          Utils.error( 'RulesParse', 'Plugin tried to change rule reference' );
      }
    });

    if( self.change( 'jsonToGroup', group, data ) != group )
      Utils.error( 'RulesParse', 'Plugin tried to change group reference' );

  }( data, this.model.root ));

  this.trigger( 'afterSetRules' );
};


/**
 * Performs value validation
 */
QueryBuilder.prototype.validateValue = function( rule, value )
{
  var validation = rule.filter.validation || {};
  var result = true;

  if( validation.callback )
    result = validation.callback.call( this, value, rule );
  else
    result = this._validateValue( rule, value );

  return this.change( 'validateValue', result, value, rule );
};

/**
 * Default validation function
 */
QueryBuilder.prototype._validateValue = function( rule, value )
{
  var filter = rule.filter;
  var operator = rule.operator;
  var validation = filter.validation || {};
  var result = true;
  var tempValue;

  if( rule.operator.nb_inputs === 1 )
    value = [value];

  for( var i = 0; i < operator.nb_inputs; i++ )
  {
    switch( filter.input )
    {
      case 'radio':
        if( value[i] === undefined || value[i].length === 0 )
        {
          if( !validation.allow_empty_value )
            result = ['radio_empty'];
        }
        break;

      case 'select':
        if( value[i] === undefined || value[i].length === 0 ||
            (filter.placeholder && value[i] == filter.placeholder_value) )
        {
          if( !validation.allow_empty_value )
            result = ['select_empty'];
        }
        break;

      default:
        tempValue = Array.isArray( value[i] ) ? value[i] : [value[i]];

        for( var j = 0; j < tempValue.length; j++ )
        {
          switch( QueryBuilder.types[filter.type] )
          {
            case 'string':
              if( tempValue[j] === undefined || tempValue[j].length === 0 )
              {
                if( !validation.allow_empty_value )
                  result = ['string_empty'];
                break;
              }
              if( validation.min !== undefined )
              {
                if( tempValue[j].length < parseInt( validation.min ) )
                {
                  result = [this.getValidationMessage( validation, 'min', 'string_exceed_min_length' ), validation.min];
                  break;
                }
              }
              if( validation.max !== undefined )
              {
                if( tempValue[j].length > parseInt( validation.max ) )
                {
                  result = [this.getValidationMessage( validation, 'max', 'string_exceed_max_length' ), validation.max];
                  break;
                }
              }
              if( validation.format )
              {
                if( typeof validation.format == 'string' )
                  validation.format = new RegExp( validation.format );

                if( !validation.format.test( tempValue[j] ) )
                {
                  result = [this.getValidationMessage( validation, 'format', 'string_invalid_format' ), validation.format];
                  break;
                }
              }
              break;

            case 'number':
              if( tempValue[j] === undefined || tempValue[j].length === 0 )
              {
                if( !validation.allow_empty_value )
                  result = ['number_nan'];
                break;
              }
              if( isNaN( tempValue[j] ) )
              {
                result = ['number_nan'];
                break;
              }
              if( filter.type == 'integer' )
              {
                if( parseInt( tempValue[j] ) != tempValue[j] )
                {
                  result = ['number_not_integer'];
                  break;
                }
              }
              else
              {
                if( parseFloat( tempValue[j] ) != tempValue[j] )
                {
                  result = ['number_not_double'];
                  break;
                }
              }
              if( validation.min !== undefined )
              {
                if( tempValue[j] < parseFloat( validation.min ) )
                {
                  result = [this.getValidationMessage( validation, 'min', 'number_exceed_min' ), validation.min];
                  break;
                }
              }
              if( validation.max !== undefined )
              {
                if( tempValue[j] > parseFloat( validation.max ) )
                {
                  result = [this.getValidationMessage( validation, 'max', 'number_exceed_max' ), validation.max];
                  break;
                }
              }
              break;
          }

          if( result !== true )
            break;
        }
    }

    if( result !== true )
      break;
  }

  return result;
};

/**
 * Returns an incremented group ID
 */
QueryBuilder.prototype.nextGroupId = function()
{
  return this.status.id + '_group_' + (this.status.group_id++);
};

/**
 * Returns an incremented rule ID
 */
QueryBuilder.prototype.nextRuleId = function()
{
  return this.status.id + '_rule_' + (this.status.rule_id++);
};

/**
 * Returns the operators for a filter
 */
QueryBuilder.prototype.getOperators = function( filter )
{
  if( typeof filter == 'string' )
    filter = this.getFilterById( filter );

  var result = [];

  for( var i = 0, l = this.operators.length; i < l; i++ )
  {
    if( filter.operators )
    {
      if( filter.operators.indexOf( this.operators[i].type ) == -1 )
        continue;
    }
    else if( this.operators[i].apply_to.indexOf( QueryBuilder.types[filter.type] ) == -1 )
    {
      continue;
    }

    result.push( this.operators[i] );
  }

  // keep sort order defined for the filter
  if( filter.operators )
  {
    result.sort( function( a, b ) {
      return filter.operators.indexOf( a.type ) - filter.operators.indexOf( b.type );
    });
  }

  return this.change( 'getOperators', result, filter );
};

/**
 * Returns a particular filter by its id
 */
QueryBuilder.prototype.getFilterById = function( id, doThrow )
{
  if( id == '-1' )
    return null;

  for( var i = 0, l = this.filters.length; i < l; i++ )
  {
    if( this.filters[i].id == id )
      return this.filters[i];
  }

  Utils.error( doThrow !== false, 'UndefinedFilter', 'Undefined filter "{0}"', id );
  return null;
};

/**
 * Returns a particular operator by its type
 */
QueryBuilder.prototype.getOperatorByType = function( type, doThrow )
{
  if( type == '-1' )
    return null;

  for( var i = 0, l = this.operators.length; i < l; i++ )
  {
    if( this.operators[i].type == type )
      return this.operators[i];
  }

  Utils.error( doThrow !== false, 'UndefinedOperator', 'Undefined operator "{0}"', type );
  return null;
};

/**
 * Returns rule's current input value
 */
QueryBuilder.prototype.getRuleInputValue = function( rule )
{
  var filter = rule.filter;
  var operator = rule.operator;
  var value = [];

  var valueContainer = rule.el.querySelector( QueryBuilder.selectors.value_container );

  for( var i = 0; i < operator.nb_inputs; i++ )
  {
    var name = Utils.escapeElementId( rule.id + '_value_' + i );

    switch( filter.input )
    {
      case 'radio':
        var checkedRadio = valueContainer.querySelector( '[name=' + name + ']:checked' );
        value.push( checkedRadio ? checkedRadio.value : undefined );
        break;

      case 'select':
        var selectEl = valueContainer.querySelector( '[name=' + name + ']' );
        if( selectEl )
        {
          var selectedOpt = selectEl.querySelector( 'option:checked' );
          value.push( selectedOpt ? selectedOpt.value : undefined );
        }
        break;

      default:
        var inputEl = valueContainer.querySelector( '[name=' + name + ']' );
        value.push( inputEl ? inputEl.value : undefined );
    }
  }

  value = value.map( function( val ) {
    if( Array.isArray( val ) )
    {
      return val.map( function( subval ) {
        return Utils.changeType( subval, filter.type );
      });
    }
    else
    {
      return Utils.changeType( val, filter.type );
    }
  });

  if( operator.nb_inputs === 1 )
    value = value[0];

  return this.change( 'getRuleValue', value, rule );
};

/**
 * Sets the value of a rule's input
 */
QueryBuilder.prototype.setRuleInputValue = function( rule, value )
{
  var filter = rule.filter;
  var operator = rule.operator;

  if( !filter || !operator )
    return;

  rule._updating_input = true;

  var valueContainer = rule.el.querySelector( QueryBuilder.selectors.value_container );

  if( operator.nb_inputs == 1 )
    value = [value];

  for( var i = 0; i < operator.nb_inputs; i++ )
  {
    var name = Utils.escapeElementId( rule.id + '_value_' + i );

    switch( filter.input )
    {
      case 'radio':
        var radioEl = valueContainer.querySelector( '[name=' + name + '][value="' + value[i] + '"]' );
        if( radioEl )
        {
          radioEl.checked = true;
          radioEl.dispatchEvent( new Event( 'change', { bubbles: true } ) );
        }
        break;

      default:
        var inputEl = valueContainer.querySelector( '[name=' + name + ']' );
        if( inputEl )
        {
          inputEl.value = value[i];
          inputEl.dispatchEvent( new Event( 'change', { bubbles: true } ) );
        }
        break;
    }
  }

  rule._updating_input = false;
};

/**
 * Parses rule flags
 */
QueryBuilder.prototype.parseRuleFlags = function( rule )
{
  var flags = extend( {}, this.settings.default_rule_flags );

  if( rule.readonly )
  {
    extend( flags, {
      filter_readonly: true,
      operator_readonly: true,
      value_readonly: true,
      no_delete: true
    });
  }

  if( rule.flags )
    extend( flags, rule.flags );

  return this.change( 'parseRuleFlags', flags, rule );
};

/**
 * Gets a copy of flags of a rule
 */
QueryBuilder.prototype.getRuleFlags = function( flags, all )
{
  if( all )
    return extend( {}, flags );

  var ret = {};
  var defaults = this.settings.default_rule_flags;
  Object.keys( defaults ).forEach( function( key ) {
    if( flags[key] !== defaults[key] )
      ret[key] = flags[key];
  });
  return ret;
};

/**
 * Parses group flags
 */
QueryBuilder.prototype.parseGroupFlags = function( group )
{
  var flags = extend( {}, this.settings.default_group_flags );

  if( group.readonly )
  {
    extend( flags, {
      condition_readonly: true,
      no_add_rule: true,
      no_add_group: true,
      no_delete: true
    });
  }

  if( group.flags )
    extend( flags, group.flags );

  return this.change( 'parseGroupFlags', flags, group );
};

/**
 * Gets a copy of flags of a group
 */
QueryBuilder.prototype.getGroupFlags = function( flags, all )
{
  if( all )
    return extend( {}, flags );

  var ret = {};
  var defaults = this.settings.default_group_flags;
  Object.keys( defaults ).forEach( function( key ) {
    if( flags[key] !== defaults[key] )
      ret[key] = flags[key];
  });
  return ret;
};

/**
 * Translate a label
 */
QueryBuilder.prototype.translate = function( category, key )
{
  if( !key )
  {
    key = category;
    category = undefined;
  }

  var translation;
  if( typeof key === 'object' )
    translation = key[this.settings.lang_code] || key['en'];
  else
    translation = (category ? this.lang[category] : this.lang)[key] || key;

  return this.change( 'translate', translation, key, category );
};

/**
 * Returns a validation message
 */
QueryBuilder.prototype.getValidationMessage = function( validation, type, def )
{
  return validation.messages && validation.messages[type] || def;
};


// =========================================================================
// Templates as plain JS functions
// =========================================================================

QueryBuilder.templates.group = function( data )
{
  var h = '<div id="' + data.group_id + '" class="rules-group-container">';
  h += '<div class="rules-group-header">';

  // Group action buttons
  h += '<div class="btn-group pull-right group-actions">';
  h += '<button type="button" class="btn btn-xs btn-success" data-add="rule">';
  h += '<i class="' + data.icons.add_rule + '"></i> ' + data.translate( "add_rule" );
  h += '</button>';

  if( data.settings.allow_groups === -1 || data.settings.allow_groups >= data.level )
  {
    h += ' <button type="button" class="btn btn-xs btn-success" data-add="group">';
    h += '<i class="' + data.icons.add_group + '"></i> ' + data.translate( "add_group" );
    h += '</button>';
  }

  if( data.level > 1 )
  {
    h += ' <button type="button" class="btn btn-xs btn-danger" data-delete="group">';
    h += '<i class="' + data.icons.remove_group + '"></i> ' + data.translate( "delete_group" );
    h += '</button>';
  }

  h += '</div>';

  // Condition radios
  h += '<div class="btn-group group-conditions">';
  data.conditions.forEach( function( condition ) {
    h += '<label class="btn btn-xs btn-primary">';
    h += '<input type="radio" name="' + data.group_id + '_cond" value="' + condition + '"> ';
    h += data.translate( "conditions", condition );
    h += '</label>';
  });
  h += '</div>';

  // Error container
  if( data.settings.display_errors )
  {
    h += '<div class="error-container"><i class="' + data.icons.error + '"></i></div>';
  }

  h += '</div>';
  h += '<div class="rules-group-body"><div class="rules-list"></div></div>';
  h += '</div>';

  return h;
};

QueryBuilder.templates.rule = function( data )
{
  var h = '<div id="' + data.rule_id + '" class="rule-container">';
  h += '<div class="rule-header">';
  h += '<div class="btn-group pull-right rule-actions">';
  h += '<button type="button" class="btn btn-xs btn-danger" data-delete="rule">';
  h += '<i class="' + data.icons.remove_rule + '"></i> ' + data.translate( "delete_rule" );
  h += '</button>';
  h += '</div>';
  h += '</div>';

  if( data.settings.display_errors )
  {
    h += '<div class="error-container"><i class="' + data.icons.error + '"></i></div>';
  }

  h += '<div class="rule-filter-container"></div>';
  h += '<div class="rule-operator-container"></div>';
  h += '<div class="rule-value-container"></div>';
  h += '</div>';

  return h;
};

QueryBuilder.templates.filterSelect = function( data )
{
  var h = '<select class="form-control" name="' + data.rule.id + '_filter">';

  if( data.settings.display_empty_filter )
    h += '<option value="-1">' + data.settings.select_placeholder + '</option>';

  data.filters.forEach( function( filter ) {
    h += '<option value="' + filter.id + '">' + data.translate( filter.label ) + '</option>';
  });

  h += '</select>';
  return h;
};

QueryBuilder.templates.operatorSelect = function( data )
{
  var h = '';

  if( data.operators.length === 1 )
  {
    h += '<span>';
    h += data.translate( "operators", data.operators[0].type );
    h += '</span>';
  }

  h += '<select class="form-control';
  if( data.operators.length === 1 ) h += ' hide';
  h += '" name="' + data.rule.id + '_operator">';

  data.operators.forEach( function( operator ) {
    h += '<option value="' + operator.type + '">';
    h += data.translate( "operators", operator.type );
    h += '</option>';
  });

  h += '</select>';
  return h;
};

QueryBuilder.templates.ruleValueSelect = function( data )
{
  var h = '<select class="form-control" name="' + data.name + '">';

  if( data.rule.filter.placeholder )
  {
    h += '<option value="' + data.rule.filter.placeholder_value + '" disabled selected>';
    h += data.rule.filter.placeholder;
    h += '</option>';
  }

  data.rule.filter.values.forEach( function( entry ) {
    h += '<option value="' + entry.value + '">' + entry.label + '</option>';
  });

  h += '</select>';
  return h;
};


/**
 * Returns group's HTML
 */
QueryBuilder.prototype.getGroupTemplate = function( group_id, level )
{
  var h = this.templates.group({
    builder: this,
    group_id: group_id,
    level: level,
    conditions: this.settings.conditions,
    icons: this.icons,
    settings: this.settings,
    translate: this.translate.bind( this )
  });

  return this.change( 'getGroupTemplate', h, level );
};

/**
 * Returns rule's HTML
 */
QueryBuilder.prototype.getRuleTemplate = function( rule_id )
{
  var h = this.templates.rule({
    builder: this,
    rule_id: rule_id,
    icons: this.icons,
    settings: this.settings,
    translate: this.translate.bind( this )
  });

  return this.change( 'getRuleTemplate', h );
};

/**
 * Returns rule's filter HTML
 */
QueryBuilder.prototype.getRuleFilterSelect = function( rule, filters )
{
  var h = this.templates.filterSelect({
    builder: this,
    rule: rule,
    filters: filters,
    icons: this.icons,
    settings: this.settings,
    translate: this.translate.bind( this )
  });

  return this.change( 'getRuleFilterSelect', h, rule, filters );
};

/**
 * Returns rule's operator HTML
 */
QueryBuilder.prototype.getRuleOperatorSelect = function( rule, operators )
{
  var h = this.templates.operatorSelect({
    builder: this,
    rule: rule,
    operators: operators,
    icons: this.icons,
    settings: this.settings,
    translate: this.translate.bind( this )
  });

  return this.change( 'getRuleOperatorSelect', h, rule, operators );
};

/**
 * Returns the rule's value select HTML
 */
QueryBuilder.prototype.getRuleValueSelect = function( name, rule )
{
  var h = this.templates.ruleValueSelect({
    builder: this,
    name: name,
    rule: rule,
    icons: this.icons,
    settings: this.settings,
    translate: this.translate.bind( this )
  });

  return this.change( 'getRuleValueSelect', h, name, rule );
};

/**
 * Returns the rule's value HTML
 */
QueryBuilder.prototype.getRuleInput = function( rule, value_id )
{
  var filter = rule.filter;
  var validation = rule.filter.validation || {};
  var name = rule.id + '_value_' + value_id;
  var c = filter.vertical ? ' class=block' : '';
  var h = '';

  switch( filter.input )
  {
    case 'radio':
      Utils.iterateOptions( filter.values, function( key, val ) {
        h += '<label' + c + '><input type="radio" name="' + name + '" value="' + key + '"> ' + val + '</label> ';
      });
      break;

    case 'select':
      h = this.getRuleValueSelect( name, rule );
      break;

    case 'number':
      h += '<input class="form-control" type="number" name="' + name + '"';
      if( validation.step !== undefined ) h += ' step="' + validation.step + '"';
      if( validation.min !== undefined ) h += ' min="' + validation.min + '"';
      if( validation.max !== undefined ) h += ' max="' + validation.max + '"';
      if( filter.placeholder ) h += ' placeholder="' + filter.placeholder + '"';
      if( filter.size ) h += ' size="' + filter.size + '"';
      h += '>';
      break;

    default: // text
      h += '<input class="form-control" type="text" name="' + name + '"';
      if( filter.placeholder ) h += ' placeholder="' + filter.placeholder + '"';
      if( filter.type === 'string' && validation.min !== undefined ) h += ' minlength="' + validation.min + '"';
      if( filter.type === 'string' && validation.max !== undefined ) h += ' maxlength="' + validation.max + '"';
      if( filter.size ) h += ' size="' + filter.size + '"';
      h += '>';
  }

  return this.change( 'getRuleInput', h, rule, name );
};


// =========================================================================
// Utils
// =========================================================================

QueryBuilder.utils = Utils;

/**
 * Iterates over radio/select options
 */
Utils.iterateOptions = function( options, tpl )
{
  if( !options )
    return;

  if( Array.isArray( options ) )
  {
    options.forEach( function( entry ) {
      if( isPlainObject( entry ) )
      {
        if( 'value' in entry )
        {
          tpl( entry.value, entry.label || entry.value );
        }
        else
        {
          var keys = Object.keys( entry );
          if( keys.length > 0 )
            tpl( keys[0], entry[keys[0]] );
        }
      }
      else
      {
        tpl( entry, entry );
      }
    });
  }
  else
  {
    Object.keys( options ).forEach( function( key ) {
      tpl( key, options[key] );
    });
  }
};

/**
 * Replaces {0}, {1}, ... in a string
 */
Utils.fmt = function( str, args )
{
  if( !Array.isArray( args ) )
    args = Array.prototype.slice.call( arguments, 1 );

  return str.replace( /{([0-9]+)}/g, function( m, i ) {
    return args[parseInt( i )];
  });
};

/**
 * Throws an Error object with custom name or logs an error
 */
Utils.error = function()
{
  var i = 0;
  var doThrow = typeof arguments[i] === 'boolean' ? arguments[i++] : true;
  var type = arguments[i++];
  var message = arguments[i++];
  var args = Array.isArray( arguments[i] ) ? arguments[i] : Array.prototype.slice.call( arguments, i );

  if( doThrow )
  {
    var err = new Error( Utils.fmt( message, args ) );
    err.name = type + 'Error';
    err.args = args;
    throw err;
  }
  else
  {
    console.error( type + 'Error: ' + Utils.fmt( message, args ) );
  }
};

/**
 * Changes the type of a value to int, float or bool
 */
Utils.changeType = function( value, type )
{
  if( value === '' || value === undefined )
    return undefined;

  switch( type )
  {
    case 'integer':
      if( typeof value === 'string' && !/^-?\d+$/.test( value ) )
        return value;
      return parseInt( value );
    case 'double':
      if( typeof value === 'string' && !/^-?\d+\.?\d*$/.test( value ) )
        return value;
      return parseFloat( value );
    default:
      return value;
  }
};

/**
 * Escapes a string for use in HTML element id selectors
 */
Utils.escapeElementId = function( str )
{
  return str ? str.replace( /(\\)?([:.\[\],])/g,
    function( $0, $1, $2 ) { return $1 ? $0 : '\\' + $2; }) : str;
};

/**
 * Expose isEmptyObject and isPlainObject on Utils for external use
 */
Utils.isEmptyObject = isEmptyObject;
Utils.isPlainObject = isPlainObject;


// =========================================================================
// Defines properties on Node prototype with getter and setter
// =========================================================================

Utils.defineModelProperties = function( obj, fields )
{
  fields.forEach( function( field ) {
    Object.defineProperty( obj.prototype, field, {
      enumerable: true,
      get: function() {
        return this.__[field];
      },
      set: function( value ) {
        var previousValue = (this.__[field] !== null && typeof this.__[field] == 'object') ?
          extend( {}, this.__[field] ) :
          this.__[field];

        this.__[field] = value;

        if( this.model !== null )
          this.model.trigger( 'update', this, field, value, previousValue );
      }
    });
  });
};


// =========================================================================
// Model
// =========================================================================

function Model()
{
  this.root = null;
  MiniEmitter.call( this );
}

Model.prototype = Object.create( MiniEmitter.prototype );
Model.prototype.constructor = Model;


// =========================================================================
// Node (base class for Group and Rule)
// =========================================================================

var Node = function( parent, el )
{
  if( !(this instanceof Node) )
    return new Node( parent, el );

  Object.defineProperty( this, '__', { value: {} } );

  el._queryBuilderModel = this;

  this.__.level = 1;
  this.__.error = null;
  this.__.flags = {};
  this.__.data = undefined;

  /**
   * DOM element (raw)
   */
  this.el = el;

  /**
   * Element ID
   */
  this.id = el.id;

  /**
   * Reference to Model
   */
  this.model = null;

  /**
   * Parent group
   */
  this.parent = parent;
};

Utils.defineModelProperties( Node, ['level', 'error', 'data', 'flags'] );

Object.defineProperty( Node.prototype, 'parent', {
  enumerable: true,
  get: function() {
    return this.__.parent;
  },
  set: function( value ) {
    this.__.parent = value;
    this.level = value === null ? 1 : value.level + 1;
    this.model = value === null ? null : value.model;
  }
});

Node.prototype.isRoot = function()
{
  return (this.level === 1);
};

Node.prototype.getPos = function()
{
  if( this.isRoot() )
    return -1;
  return this.parent.getNodePos( this );
};

Node.prototype.drop = function()
{
  var model = this.model;

  if( !!this.parent )
    this.parent.removeNode( this );

  delete this.el._queryBuilderModel;

  if( model !== null )
    model.trigger( 'drop', this );
};

Node.prototype.moveAfter = function( target )
{
  if( !this.isRoot() )
    this.move( target.parent, target.getPos() + 1 );
};

Node.prototype.moveAtBegin = function( target )
{
  if( !this.isRoot() )
  {
    if( target === undefined )
      target = this.parent;
    this.move( target, 0 );
  }
};

Node.prototype.moveAtEnd = function( target )
{
  if( !this.isRoot() )
  {
    if( target === undefined )
      target = this.parent;
    this.move( target, target.length() === 0 ? 0 : target.length() - 1 );
  }
};

Node.prototype.move = function( target, index )
{
  if( !this.isRoot() )
  {
    if( typeof target === 'number' )
    {
      index = target;
      target = this.parent;
    }

    this.parent.removeNode( this );
    target.insertNode( this, index, false );

    if( this.model !== null )
      this.model.trigger( 'move', this, target, index );
  }
};


// =========================================================================
// Group
// =========================================================================

var Group = function( parent, el )
{
  if( !(this instanceof Group) )
    return new Group( parent, el );

  Node.call( this, parent, el );

  this.rules = [];
  this.__.condition = null;
};

Group.prototype = Object.create( Node.prototype );
Group.prototype.constructor = Group;

Utils.defineModelProperties( Group, ['condition'] );

Group.prototype.empty = function()
{
  this.each( 'reverse', function( rule ) {
    rule.drop();
  }, function( group ) {
    group.drop();
  });
};

Group.prototype.drop = function()
{
  this.empty();
  Node.prototype.drop.call( this );
};

Group.prototype.length = function()
{
  return this.rules.length;
};

Group.prototype.insertNode = function( node, index, trigger )
{
  if( index === undefined )
    index = this.length();

  this.rules.splice( index, 0, node );
  node.parent = this;

  if( trigger && this.model !== null )
    this.model.trigger( 'add', this, node, index );

  return node;
};

Group.prototype.addGroup = function( el, index )
{
  return this.insertNode( new Group( this, el ), index, true );
};

Group.prototype.addRule = function( el, index )
{
  return this.insertNode( new Rule( this, el ), index, true );
};

Group.prototype.removeNode = function( node )
{
  var index = this.getNodePos( node );
  if( index !== -1 )
  {
    node.parent = null;
    this.rules.splice( index, 1 );
  }
};

Group.prototype.getNodePos = function( node )
{
  return this.rules.indexOf( node );
};

Group.prototype.each = function( reverse, cbRule, cbGroup, context )
{
  if( typeof reverse !== 'boolean' && typeof reverse !== 'string' )
  {
    context = cbGroup;
    cbGroup = cbRule;
    cbRule = reverse;
    reverse = false;
  }
  context = context === undefined ? null : context;

  var i = reverse ? this.rules.length - 1 : 0;
  var l = reverse ? 0 : this.rules.length - 1;
  var c = reverse ? -1 : 1;
  var next = function() {
    return reverse ? i >= l : i <= l;
  };
  var stop = false;

  for( ; next(); i += c )
  {
    if( this.rules[i] instanceof Group )
    {
      if( !!cbGroup )
        stop = cbGroup.call( context, this.rules[i] ) === false;
    }
    else if( !!cbRule )
    {
      stop = cbRule.call( context, this.rules[i] ) === false;
    }

    if( stop )
      break;
  }

  return !stop;
};

Group.prototype.contains = function( node, recursive )
{
  if( this.getNodePos( node ) !== -1 )
    return true;
  else if( !recursive )
    return false;
  else
  {
    return !this.each( function() {
      return true;
    }, function( group ) {
      return !group.contains( node, true );
    });
  }
};


// =========================================================================
// Rule
// =========================================================================

var Rule = function( parent, el )
{
  if( !(this instanceof Rule) )
    return new Rule( parent, el );

  Node.call( this, parent, el );

  this._updating_value = false;
  this._updating_input = false;

  this.__.filter = null;
  this.__.operator = null;
  this.__.value = undefined;
};

Rule.prototype = Object.create( Node.prototype );
Rule.prototype.constructor = Rule;

Utils.defineModelProperties( Rule, ['filter', 'operator', 'value'] );

Rule.prototype.isRoot = function()
{
  return false;
};


// =========================================================================
// Expose classes
// =========================================================================

QueryBuilder.Group = Group;
QueryBuilder.Rule = Rule;


// =========================================================================
// English translations (inline)
// =========================================================================

QueryBuilder.regional['en'] = {
  "__locale": "English (en)",
  "add_rule": "Add rule",
  "add_group": "Add group",
  "delete_rule": "Delete",
  "delete_group": "Delete",
  "conditions": {
    "AND": "AND",
    "OR": "OR"
  },
  "operators": {
    "equal": "equal",
    "not_equal": "not equal",
    "in": "in",
    "not_in": "not in",
    "less": "less",
    "less_or_equal": "less or equal",
    "greater": "greater",
    "greater_or_equal": "greater or equal",
    "between": "between",
    "not_between": "not between",
    "begins_with": "begins with",
    "not_begins_with": "doesn't begin with",
    "contains": "contains",
    "not_contains": "doesn't contain",
    "ends_with": "ends with",
    "not_ends_with": "doesn't end with",
    "is_empty": "is empty",
    "is_not_empty": "is not empty",
    "is_null": "is null",
    "is_not_null": "is not null"
  },
  "errors": {
    "no_filter": "No filter selected",
    "empty_group": "The group is empty",
    "radio_empty": "No value selected",
    "checkbox_empty": "No value selected",
    "select_empty": "No value selected",
    "string_empty": "Empty value",
    "string_exceed_min_length": "Must contain at least {0} characters",
    "string_exceed_max_length": "Must not contain more than {0} characters",
    "string_invalid_format": "Invalid format ({0})",
    "number_nan": "Not a number",
    "number_not_integer": "Not an integer",
    "number_not_double": "Not a real number",
    "number_exceed_min": "Must be greater than {0}",
    "number_exceed_max": "Must be lower than {0}",
    "number_wrong_step": "Must be a multiple of {0}",
    "number_between_invalid": "Invalid values, {0} is greater than {1}",
    "boolean_not_valid": "Not a boolean",
    "operator_not_multiple": "Operator \"{1}\" cannot accept multiple values"
  }
};

QueryBuilder.defaults({ lang_code: 'en' });


// =========================================================================
// Factory function
// =========================================================================

/**
 * Creates a QueryBuilder instance on a given element.
 * @param {HTMLElement|string} element - DOM element or element ID
 * @param {object} options
 * @returns {QueryBuilder}
 */
QueryBuilder.create = function( element, options )
{
  if( typeof element === 'string' )
    element = document.getElementById( element );

  if( !element )
    Utils.error( 'Config', 'No target element found' );

  options = options || {};

  // Bind events from options before init
  var events = options.events;
  delete options.events;

  var builder = new QueryBuilder( element, options );

  // Store instance on element for external access
  element._queryBuilder = builder;

  // Bind events before init so they catch initialization events
  if( events )
  {
    Object.keys( events ).forEach( function( eventName ) {
      // Support space-separated event names
      eventName.split( ' ' ).forEach( function( name ) {
        if( name.trim() )
          builder.on( name.trim(), events[name.trim()] || events[eventName] );
      });
    });
  }

  builder.init( options.rules );

  return builder;
};


// =========================================================================
// Expose globally
// =========================================================================

root.QueryBuilder = QueryBuilder;

}( this ));

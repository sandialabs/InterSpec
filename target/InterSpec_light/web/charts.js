/* InterSpec Light - the spectrum and time-series charts, and routing their events to the app. */
'use strict';

App.initCharts = function() {
  const spec = new SpectrumChartD3( 'spectrum-chart', {
    xlabel: 'Energy (keV)',
    ylabel: 'Counts',
    yscale: 'log',
    showLegend: true,
    showMouseStats: true,
    showRefLineInfoForMouseOver: true,
    allowDragRoiExtent: true,
    compactXAxis: false,
    showXAxisSliderChart: false,
    refLineVerbosity: 1
  } );
  App.spectrumChart = spec;

  const time = new D3TimeChart( 'time-chart', {
    xtitle: 'Time (s)',
    y1title: 'Gamma CPS',
    y2title: 'Neutron CPS',
    compactXAxis: true,
    chartLineWidth: 1.0
  } );
  App.timeChart = time;

  // The charts do not watch their own size
  const observer = new ResizeObserver( (entries) => {
    for( const entry of entries )
    {
      if( entry.target.id === 'spectrum-chart' )
        spec.handleResize();
      else if( (entry.target.id === 'time-chart') && !entry.target.hidden )
        time.handleResize();
    }
  } );
  observer.observe( document.getElementById( 'spectrum-chart' ) );
  observer.observe( document.getElementById( 'time-chart' ) );

  // Y-axis type and energy strip chart; the chart can also toggle these itself (clicking the
  //  axis titles), so it reports changes back, and the choice is remembered per browser.
  const yscale = document.getElementById( 'yscale' );
  const slider = document.getElementById( 'show-slider' );
  const remember = ( key, value ) => { try{ localStorage.setItem( key, value ); }catch( e ){} };
  const recall = ( key ) => { try{ return localStorage.getItem( key ); }catch( e ){ return null; } };

  yscale.addEventListener( 'change', () => {
    spec.setYAxisType( yscale.value );
    remember( 'interspec-light-yscale', yscale.value );
  } );
  slider.addEventListener( 'change', () => {
    spec.setShowXAxisSliderChart( slider.checked );
    remember( 'interspec-light-slider', slider.checked ? '1' : '0' );
  } );

  const savedY = recall( 'interspec-light-yscale' );
  if( (savedY === 'lin') || (savedY === 'log') || (savedY === 'sqrt') )
    yscale.value = savedY;
  spec.setYAxisType( yscale.value );
  slider.checked = (recall( 'interspec-light-slider' ) === '1');
  spec.setShowXAxisSliderChart( slider.checked );

  // Replies to the chart are deferred, so we never re-enter it from inside its own event handler.
  const later = (fn) => setTimeout( fn, 0 );

  App.chartHandlers['spectrum-chart'] = {
    yAxisTypeChanged: (ytype) => {
      yscale.value = ytype;
      remember( 'interspec-light-yscale', ytype );
    },
    sliderChartDisplayed: (show) => {
      slider.checked = !!show;
      remember( 'interspec-light-slider', show ? '1' : '0' );
    },
    doubleclicked: (energy, count, refLineParent) => {
      App.peaks.fitAt( energy, refLineParent );
    },
    rightclicked: (energy, count, pageX, pageY, refLineParent) => {
      App.peaks.showContextMenu( energy, pageX, pageY );
    },
    shiftkeydragged: (e0, e1) => {
      App.run( 'erasePeaks', { e0: e0, e1: e1 }, { busy: true } );
    },
    roiDrag: (newLower, newUpper, px, origLower, specType, isFinal) => {
      if( specType !== 'FOREGROUND' )
      {
        later( () => spec.updateRoiBeingDragged( null ) );
        return;
      }
      let result = null;
      try
      {
        result = App.call( 'roiDrag', { newLower: newLower, newUpper: newUpper, px: px,
                                         origLower: origLower, isFinal: !!isFinal } );
      }catch( e )
      {
        App.toast( e.message, true );
      }
      later( () => {
        if( !isFinal )
          spec.updateRoiBeingDragged( result ? result.roiDragPeaks : null );
        else
          App.applyState( result );
      } );
    },
    fitRoiDrag: (lower, upper, npeaks, isFinal) => {
      let result = null;
      try
      {
        result = App.call( 'fitRoiDrag', { lower: lower, upper: upper, npeaks: npeaks, isFinal: !!isFinal } );
      }catch( e )
      {
        App.toast( e.message, true );
      }
      later( () => {
        if( !isFinal )
          spec.updateRoiBeingDragged( result ? result.roiDragPeaks : null );
        else
          App.applyState( result );
      } );
    }
  };

  App.chartHandlers['time-chart'] = {
    timedragged: (first, last, mods) => {
      App.run( 'timeDrag', { first: first, last: last, mods: mods }, { busy: true } );
    }
  };

  App.renderers.push( App.renderCharts );
};


/** Pixels per keV of the spectrum chart, for the peak-fit heuristics. */
App.pixelsPerKeV = function() {
  const chart = App.spectrumChart;
  if( !chart || !chart.xScale )
    return 2.0;
  const d = chart.xScale.domain(), r = chart.xScale.range();
  const ppk = (r[1] - r[0]) / (d[1] - d[0]);
  return (isFinite( ppk ) && (ppk > 0)) ? ppk : 2.0;
};


App.renderCharts = function( state ) {
  const spec = App.spectrumChart;

  if( state.spectra )
  {
    const hasFore = state.spectra.some( (s) => s.type === 'FOREGROUND' );
    App.hasForeground = hasFore;
    document.getElementById( 'empty-msg' ).hidden = hasFore;
    spec.setData( state.spectra.length ? { spectra: state.spectra } : null, !!state.resetDomain );
  }else if( state.peaks )
  {
    spec.setRoiData( state.peaks, 'FOREGROUND' );
  }

  if( 'timeChart' in state )
  {
    const div = document.getElementById( 'time-chart' );
    const show = !!state.timeChart;
    const wasHidden = div.hidden;
    div.hidden = !show;
    if( show )
    {
      App.timeChart.setData( state.timeChart );
      App.timeChart.setHighlightRegions( (state.timeHighlights && state.timeHighlights.length) ? state.timeHighlights : null );
      if( wasHidden )
        requestAnimationFrame( () => App.timeChart.handleResize() );
    }
    if( wasHidden === show )
      requestAnimationFrame( () => spec.handleResize() );
  }
};

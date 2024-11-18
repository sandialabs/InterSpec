# Leaflet.Measure 
Leaflet.Measure Is a leaflet plugin for measuring distances and areas.  

Online [DEMO](https://ptma.github.io/Leaflet.Measure/examples/measure.html).

## Example
```javascript
// 1. add control to map
var map = L.map("map", {
        center: [29, 120],
    });
L.control.measure().addTo(map);

// 2. using action directly
var measureAction = new L.MeasureAction(map, {
    model: "distance", // 'area' or 'distance', default is 'distance'
});
// measureAction.setModel('area');
measureAction.enable();
```

## L.Control.Measure 
### Options
| Option | Type | Default | Description |
|--------|------|---------|-------------|
| position | String | 'topleft' | The position of the control. |
| title | String | 'Measurement' | The title of the control. |
| collapsed | Boolean | false | If true, the control will be collapsed into an icon and expanded on mouse hover or touch. |
| color | String | '#FF0080'| The color of the lines or polygons. |

## L.MeasureAction 
### Options
| Option | Type | Default | Description |
|--------|------|---------|-------------|
| model | String | 'distance' | Measurement mode. Possible values are 'distance' or 'area'. |
| color | String | '#FF0080'| The color of the lines or polygons. |

### Methods
| Method | Returns | Description |
|--------|---------|-------------|
| setModel(\<String\> model) | this | Sets the measurement mode. Possible values are 'distance' or 'area'. |

## Customize Language
```javascript
L.Measure = {
    linearMeasurement: "Distance measurement",
    areaMeasurement: "Area measurement",
    start: "Start",
    meter: "m",
    meterDecimals: 0,
    kilometer: "km",
    kilometerDecimals: 2,
    squareMeter: "m²",
    squareMeterDecimals: 0,
    squareKilometers: "km²",
    squareKilometersDecimals: 2
};
```

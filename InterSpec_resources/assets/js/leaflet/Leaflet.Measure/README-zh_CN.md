# Leaflet.Measure 
Leaflet.Measure 是一个 Leaflet 的测量距离和面积的工具插件。  

在线 [示例](https://ptma.gitee.io/leaflet.measure/examples/measure.html).

## 用法示例
```javascript
// 1. 添加控件
var map = L.map("map", {
        center: [29, 120],
    });
    
L.control.measure().addTo(map);

// 2. 直接执行动作
var measureAction = new L.MeasureAction(map, {
    model: "distance", // 'area' or 'distance', default is 'distance'
});
// measureAction.setModel('area');
measureAction.enable();
```

## 控件（L.Control.Measure）
### 选项
| 选项 | 类型 | 默认值 | 描述 |
|--------|------|---------|-------------|
| position | String | 'topleft' | 控件位置。 |
| title | String | 'Legend' | 控件面板的标题。 |
| collapsed | Boolean | false | 面板是否默认展开。 |
| color | String | '#FF0080'| 测量线条的颜色。 |

## 动作（L.MeasureAction）
### 选项
| 选项 | 类型 | 默认值 | 描述 |
|--------|------|---------|-------------|
| model | String | 'distance' | 测量模式. 可以是'distance' 或 'area'. |
| color | String | '#FF0080'| The color of the lines or polygons. |

### 方法
| 方法 | 返回 | 描述 |
|--------|---------|-------------|
| setModel(\<String\> model) | this | 设置测量模式. 可以是'distance' 或 'area'. |

## 显示语言定义
```javascript
L.Measure = {
    linearMeasurement: "距离测量",
    areaMeasurement: "面积测量",
    start: "开始",
    meter: "米",
    meterDecimals: 0,
    kilometer: "公里",
    kilometerDecimals: 2,
    squareMeter: "平方米",
    squareMeterDecimals: 0,
    squareKilometers: "平方公里",
    squareKilometersDecimals: 2
};
```

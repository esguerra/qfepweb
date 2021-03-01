

var DIR = "imgs/";

var nodes = null;
var edges = null;
var network = null;

// Called when the Visualization API is loaded.
function draw() {
  // create people.
  // value corresponds with the age of the person
  var DIR = "imgs/";
//  nodes = [
//    { id: 1, shape: "circularImage", image: DIR + "17.png" },
//    { id: 2, shape: "circularImage", image: DIR + "1h1q.png" },
//    { id: 3, shape: "circularImage", image: DIR + "1h1r.png" },
//    { id: 5, shape: "circularImage", image: DIR + "1h1s.png" },
//    { id: 6, shape: "circularImage", image: DIR + "1oiu.png" },
//    { id: 7, shape: "circularImage", image: DIR + "1oiy.png" },
//    { id: 8, shape: "circularImage", image: DIR + "20.png" },
//    { id: 9, shape: "circularImage", image: DIR + "21.png" },
//    { id: 10, shape: "circularImage", image: DIR + "22.png" },
//    { id: 11, shape: "circularImage", image: DIR + "26.png" },
//    { id: 12, shape: "circularImage", image: DIR + "28.png" },
//    { id: 13, shape: "circularImage", image: DIR + "29.png" },
//  ];
//
//
//  // create connections between people
//  // value corresponds with the amount of contact between two people
//  edges = [
//    { from: 1, to: 2 },
//    { from: 2, to: 3 },
//    { from: 3, to: 4 },
//    { from: 4, to: 5 },
//    { from: 5, to: 6 },
//    { from: 6, to: 7 },
//    { from: 7, to: 8 },
//    { from: 8, to: 9 },
//    { from: 9, to: 10 },
//    { from: 10, to: 11 },
//    { from: 11, to: 12 },
//    { from: 12, to: 13 },
//    { from: 3, to: 6 },
//    { from: 4, to: 8 },
//    { from: 2, to: 11 },
//  ];
    
  // load input graph file.
  var json = $.getJSON("input/graph.json")
  .done(function(data){
    var data = {
      nodes: data.nodes,
      edges: data.edges
    };

  });
  // create a network
  var container = document.getElementById("mynetwork");
  var data = {
    nodes: nodes,
    edges: edges,
  };
  var options = {
    nodes: {
      // could also make write PNGs as circular images, then set shape : image and omit size: n.
      // probably should do this. node clicking is inconsistent currently.
      //shape: "image",
      borderWidth: 1,
      size: 55,
      color: {
        border: "#222222",
        background: "white",
        highlight: {
            border: '#2B7CE9',
            background: "white"
            }
      },
      font: { color: "#eeeeee" },
    },
    edges: {
      color: "lightgray",
    },
  };
            
  network = new vis.Network(container, data, options);
}

window.addEventListener("load", () => {
  draw();
});


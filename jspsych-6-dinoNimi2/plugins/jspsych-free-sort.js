/**
 * jspsych-free-sort
 * plugin for drag-and-drop sorting of a collection of images
 * Josh de Leeuw
 *
 * documentation: docs.jspsych.org
 */


jsPsych.plugins['free-sort'] = (function() {

  var plugin = {};

  jsPsych.pluginAPI.registerPreload('free-sort', 'stimuli', 'image');

  plugin.info = {
    name: 'free-sort',
    description: '',
    parameters: {
      stimuli: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Stimuli',
        default: undefined,
        array: true,
        description: 'Images to be displayed.'
      },
      stim_height: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Stimulus height',
        default: 70,
        description: 'Height of images in pixels.'
      },
      stim_width: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Stimulus width',
        default: 110,
        description: 'Width of images in pixels'
      },
      sort_area_height: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Sort area height',
        default: 800,
        description: 'The height of the container that subjects can move the stimuli in.'
      },
      sort_area_width: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Sort area width',
        default: 800,
        description: 'The width of the container that subjects can move the stimuli in.'
      },
      prompt: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Prompt',
        default: null,
        description: 'It can be used to provide a reminder about the action the subject is supposed to take.'
      },
      prompt_location: {
        type: jsPsych.plugins.parameterType.SELECT,
        pretty_name: 'Prompt location',
        options: ['above','below'],
        default: 'above',
        description: 'Indicates whether to show prompt "above" or "below" the sorting area.'
      },
      button_label: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Button label',
        default:  'Next',
        description: 'The text that appears on the button to continue to the next trial.'
      }
    }
  }

  plugin.trial = function(display_element, trial) {

    var start_time = performance.now();

    var html = "";
    // check if there is a prompt and if it is shown above
    if (trial.prompt !== null && trial.prompt_location == "above") {
      html += trial.prompt;
    }

    html += '<div '+
      'id="jspsych-free-sort-arena" '+
      'class="jspsych-free-sort-arena" '+
      'style="position: relative; width:'+trial.sort_area_width+'px; height:'+trial.sort_area_height+'px; border:2px solid #444;"'+
      '></div>';
    // check if prompt exists and if it is shown below
    if (trial.prompt !== null && trial.prompt_location == "below") {
      html += trial.prompt;
    }

    display_element.innerHTML = html;

    // draw four quadrants
    for(var qi=0;qi<4;qi++){
      //if (qi==0) {topqi=20;leftqi=20}else{topqi=(20+260+20+10);leftqi=(20+360+20+10)}
      if (qi==0 || qi==2) {leftqi=20}
      if (qi==1 || qi==3) {leftqi=20+360+20+80} //20+360+20+10
      if (qi==0 || qi==1) {topqi=20}
      if (qi==2 || qi==3) {topqi=20+260+20+80}
      if (qi==0){qiText="Family 1"}
      if (qi==1){qiText="Family 2"}
      if (qi==2){qiText="Family 3"}
      if (qi==3){qiText="Family 4"}
        // 370 270
      display_element.querySelector("#jspsych-free-sort-arena").innerHTML += '<rect '+
      'style="position: absolute; cursor: move; width:'+300+'px; height:'+200+'px; '+
      'top:'+topqi+'px; left:'+leftqi+'px; background-color: #E8E8E8; ">'+
      '<rect '+
        'style="position: absolute; cursor: move; width:'+(trial.stim_width+10)+'px; '+
        'height:'+(trial.stim_height+10)+'px; top:'+70+'px; left:'+15+'px; '+
        'background-color: white; border-style:dashed;stroke-width:10">'+
      '</rect>'+
      '<div '+
        'style="position: absolute; cursor: move; width:'+(trial.stim_width+10)+'px; '+
        'height:'+(trial.stim_height+10)+'px; top:'+70+'px; left:'+(15+trial.stim_width+30)+'px; '+
        'background-color: white; border-style:dashed;stroke-width:10">'+
      '</div>'+

      '<text font-family="Verdana" font-size="50" fill="blue"><p>'+qiText+'</text>'+
      '</rect>';
    };

    // store initial location data
    var init_locations = [];

    for (var i = 0; i < trial.stimuli.length; i++) {
      //var coords = random_coordinate(trial.sort_area_width - trial.stim_width, trial.sort_area_height - trial.stim_height);
      
      if(i==0){xcoord=20}else{xcoord=(20+(trial.stim_width/1.2)*i)}
      y_midpoint = (trial.sort_area_height/2)-(trial.stim_height/2)
      if(i%2==0){ycoord=y_midpoint-20}else{ycoord=y_midpoint+20}
  
      var coords = {x:xcoord,y:ycoord}

      display_element.querySelector("#jspsych-free-sort-arena").innerHTML += '<img '+
        'src="'+trial.stimuli[i]+'" '+
        'data-src="'+trial.stimuli[i]+'" '+
        'class="jspsych-free-sort-draggable" '+
        'draggable="false" '+
        'style="position: absolute; cursor: move; width:'+trial.stim_width+'px; height:'+trial.stim_height+'px; top:'+coords.y+'px; left:'+coords.x+'px;">'+
        '</img>';

      init_locations.push({
        "src": trial.stimuli[i],
        "x": coords.x,
        "y": coords.y
      });
    }

    display_element.innerHTML += '<p><button id="jspsych-free-sort-done-btn" class="jspsych-btn">'+trial.button_label+'</button>';

    var maxz = 1;

    var moves = [];

    var draggables = display_element.querySelectorAll('.jspsych-free-sort-draggable');

    for(var i=0;i<draggables.length; i++){
      draggables[i].addEventListener('mousedown', function(event){
        var x = event.pageX - event.currentTarget.offsetLeft;
        var y = event.pageY - event.currentTarget.offsetTop - window.scrollY;
        var elem = event.currentTarget;
        elem.style.zIndex = ++maxz;

        var mousemoveevent = function(e){
          elem.style.top =  Math.min(trial.sort_area_height - trial.stim_height, Math.max(0,(e.clientY - y))) + 'px';
          elem.style.left = Math.min(trial.sort_area_width  - trial.stim_width,  Math.max(0,(e.clientX - x))) + 'px';
        }
        document.addEventListener('mousemove', mousemoveevent);

        var mouseupevent = function(e){
          document.removeEventListener('mousemove', mousemoveevent);
          moves.push({
            "src": elem.dataset.src,
            "x": elem.offsetLeft,
            "y": elem.offsetTop
          });
          document.removeEventListener('mouseup', mouseupevent);
        }
        document.addEventListener('mouseup', mouseupevent);
      });
    }

    display_element.querySelector('#jspsych-free-sort-done-btn').addEventListener('click', function(){

      var end_time = performance.now();
      var rt = end_time - start_time;
      // gather data
      // get final position of all objects
      var final_locations = [];
      var matches = display_element.querySelectorAll('.jspsych-free-sort-draggable');
      for(var i=0; i<matches.length; i++){
        final_locations.push({
          "src": matches[i].dataset.src,
          "x": parseInt(matches[i].style.left),
          "y": parseInt(matches[i].style.top)
        });
      }
      //console.log([final_locations[1],final_locations[1].x]);

      // Have all animals been moved?
      var allMoved='yes';
      for(ax=0;ax<Object.keys(tax_bottomUp).length;ax++){
        if((final_locations[ax].x-init_locations[ax].x)==0 && (final_locations[ax].y-init_locations[ax].y)==0){
          allMoved='no'}
      };


      allElements = ["Barilium", "Charachodontosaurus", "Uteodon", 
                      "Argentinosaurus", "Saturnalia", "Isasicursor", "Velociraptor", "Macrogryphosaurus"];
      pairs = [];
      for (var i = 0; i < 4; i++) {
        Element_i = allElements[0];
        allElements.splice(0,1);
        idx_Element_i = sub_familyName_freeSort.indexOf(Element_i);
        Element_i_coords = [final_locations[idx_Element_i].x,final_locations[idx_Element_i].y];
        ij_distances = [];
        for(j=0;j<allElements.length;j++){
          Element_j = allElements[j];
          idx_Element_j = sub_familyName_freeSort.indexOf(Element_j);
          Element_j_coords = [final_locations[idx_Element_j].x,final_locations[idx_Element_j].y];
          dist_x = Math.pow(Element_i_coords[0]-Element_j_coords[0],2);
          dist_y = Math.pow(Element_i_coords[1]-Element_j_coords[1],2);
          pair_ij_euclid = Math.sqrt(dist_x + dist_y);
          ij_distances.push(pair_ij_euclid);
        }
        minDistance = Math.min.apply(Math,ij_distances);
        idx_minDistance = ij_distances.indexOf(minDistance);
        element_with_minDist = allElements[idx_minDistance];
        allElements.splice(idx_minDistance,1);
        pairs.push([Element_i,element_with_minDist]);
      }
      // check if which of the pairings are correct
      nCorrectPairs = 0;
      correctPairs_freeSort = [];
      correctFamilies_freeSort = [];
      for(i=0;i<4;i++){
        pair_i = pairs[[i]];
        pair_i_element1 = pair_i[0];
        pair_i_element2 = pair_i[1];
        pair_i_element1_familiy = tax_bottomUp[pair_i_element1];
        pair_i_element2_familiy = tax_bottomUp[pair_i_element2];
        if(pair_i_element1_familiy==pair_i_element2_familiy){
          correctPairs_freeSort.push([pair_i_element1,pair_i_element2])
          nCorrectPairs += 1;
          correctFamilies_freeSort.push(pair_i_element1_familiy);
          nCorrect_freeSort[pair_i_element1_familiy][count_freeSort] += 1
          //nCorrect_freeSort[pair_i_element2_familiy][count_freeSort] += 1
        }else{
          //nCorrect_freeSort[pair_i_element1_familiy][count_freeSort].push(0)
          //nCorrect_freeSort[pair_i_element2_familiy][count_freeSort].push(0)
        }
      }

      /*nCorrectPairs = 0;
      for(Fam_x=0;Fam_x<Object.keys(dino_pics).length;Fam_x++){
        //var element = 'Theropoda';
        var element = Object.keys(dino_pics)[Fam_x];
        idx_element = [];
        for(i=0;i<basic_familyName_freeSort.length;i++){
          if(basic_familyName_freeSort[i]==element){idx_element.push(i)}
        }
        element1_coords = [final_locations[idx_element[0]].x,final_locations[idx_element[0]].y];
        element2_coords = [final_locations[idx_element[1]].x,final_locations[idx_element[1]].y];
        dist_x = Math.pow(element1_coords[0]-element2_coords[0],2);
        dist_y = Math.pow(element1_coords[1]-element2_coords[1],2);
        pair_px_euclid = Math.sqrt(dist_x + dist_y);

        idx_remaining = [];
        for(i=0;i<basic_familyName_freeSort.length;i++){
          if(basic_familyName_freeSort[i]!=element){idx_remaining.push(i)}
        }

        element1_distances = [];
        for(i=0;i<idx_remaining.length;i++){
          ext_element_i_coords = [final_locations[idx_remaining[i]].x,final_locations[idx_remaining[i]].y];
          dist_x = Math.pow(element1_coords[0]-ext_element_i_coords[0],2);
          dist_y = Math.pow(element1_coords[1]-ext_element_i_coords[1],2);
          element1_distances.push(Math.sqrt(dist_x + dist_y));
        }
        element2_distances = [];
        for(i=0;i<idx_remaining.length;i++){
          ext_element_i_coords = [final_locations[idx_remaining[i]].x,final_locations[idx_remaining[i]].y];
          dist_x = Math.pow(element2_coords[0]-ext_element_i_coords[0],2);
          dist_y = Math.pow(element2_coords[1]-ext_element_i_coords[1],2);
          element2_distances.push(Math.sqrt(dist_x + dist_y));
        }

        // Is within-distance smaller than between-distance?
        all_between_distances = element1_distances.concat(element2_distances);
        nr_largerBetweenDistances = 0;
        for (var i = 0; i < all_between_distances.length; i++) {
          if(pair_px_euclid > all_between_distances[i]){nr_largerBetweenDistances += 1}
        }
   
        if(nr_largerBetweenDistances==0){
          nCorrectPairs += 1
          nCorrect_freeSort[element][count_freeSort].push(1)
        }else{
          nCorrect_freeSort[element][count_freeSort].push(0)
        }
      }
      console.log(nCorrectPairs);*/

      var trial_data = {
        "init_locations": JSON.stringify(init_locations),
        "moves": JSON.stringify(moves),
        "final_locations": JSON.stringify(final_locations),
        "Pairings":JSON.stringify(pairs),
        "CorrectPairings":JSON.stringify(correctPairs_freeSort),
        "CorrectFamilies":JSON.stringify(correctFamilies_freeSort),
        "nCorrectPairs": nCorrectPairs,
        "rt": rt
      };
      //allMoved='yes';
      if(allMoved=='no'){alert("Jätkamiseks pead iga looma asetama ühte kaheksast nelinurgast.");}else{
        // advance to next part
        display_element.innerHTML = '';
        jsPsych.finishTrial(trial_data);
      }
      
    });

  };

  // helper functions

  function random_coordinate(max_width, max_height) {
    var rnd_x = Math.floor(Math.random() * (max_width - 1));
    var rnd_y = Math.floor(Math.random() * (max_height - 1));

    return {
      x: rnd_x,
      y: rnd_y
    };
  }

  return plugin;
})();

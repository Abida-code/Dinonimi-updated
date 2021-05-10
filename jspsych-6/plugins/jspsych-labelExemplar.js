/**
 * jspsych-image-button-response
 * Josh de Leeuw
 *
 * plugin for displaying a stimulus and getting a keyboard response
 *
 * documentation: docs.jspsych.org
 *
 **/

jsPsych.plugins["labelExemplar"] = (function() {

  var plugin = {};

  jsPsych.pluginAPI.registerPreload('image-button-response', 'stimulus', 'image');

  plugin.info = {
    name: 'image-button-response',
    description: '',
    parameters: {
      stimulus: {
        type: jsPsych.plugins.parameterType.IMAGE,
        pretty_name: 'Stimulus',
        default: undefined,
        description: 'The image to be displayed'
      },
      species_name: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Species_name',
        default: undefined,
        description: 'The name of the selected species.'
      },
      stimulus_height: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Image height',
        default: null,
        description: 'Set the image height in pixels'
      },
      stimulus_width: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Image width',
        default: null,
        description: 'Set the image width in pixels'
      },
      maintain_aspect_ratio: {
        type: jsPsych.plugins.parameterType.BOOL,
        pretty_name: 'Maintain aspect ratio',
        default: true,
        description: 'Maintain the aspect ratio after setting width or height'
      },
      choices: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Choices',
        default: undefined,
        array: true,
        description: 'The labels for the buttons.'
      },
      button_html: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Button HTML',
        default: '<button class="jspsych-btn">%choice%</button>',
        array: true,
        description: 'The html of the button. Can create own style.'
      },
      prompt: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Prompt',
        default: null,
        description: 'Any content here will be displayed under the button.'
      },
      stimulus_duration: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Stimulus duration',
        default: null,
        description: 'How long to hide the stimulus.'
      },
      trial_duration: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Trial duration',
        default: null,
        description: 'How long to show the trial.'
      },
      margin_vertical: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Margin vertical',
        default: '0px',
        description: 'The vertical margin of the button.'
      },
      margin_horizontal: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Margin horizontal',
        default: '8px',
        description: 'The horizontal margin of the button.'
      },
      response_ends_trial: {
        type: jsPsych.plugins.parameterType.BOOL,
        pretty_name: 'Response ends trial',
        default: true,
        description: 'If true, then trial will end when user responds.'
      },
    }
  }

  plugin.trial = function(display_element, trial) {

    // display stimulus
    var html = '<img src="'+trial.stimulus+'" id="jspsych-image-button-response-stimulus" style="';
    if(trial.stimulus_height !== null){
      html += 'height:'+trial.stimulus_height+'px; '
      if(trial.stimulus_width == null && trial.maintain_aspect_ratio){
        html += 'width: auto; ';
      }
    }
    if(trial.stimulus_width !== null){
      html += 'width:'+trial.stimulus_width+'px; '
      if(trial.stimulus_height == null && trial.maintain_aspect_ratio){
        html += 'height: auto; ';
      }
    }
    html +='"></img>';
    // display species_properties
    var props = '<p><b> Name: </b>'+species_clicked+'</p>';
    //html += props;
    var PropTab = '<table id="properties" align="center"> ';
      PropTab += '<tr> <th>Dimension</th><th>Feature</th><th>Illustration</th></tr>';
      PropTab += '<tr><td><b>Size</b></td><td>'+properties_species_clicked[0]+'</td><td><img src="'+props_illustrated_species_clicked[0]+'"></td></tr>';
      PropTab += '<tr><td><b>Hip</b></td><td>'+properties_species_clicked[1]+'</td><td><img src="'+props_illustrated_species_clicked[1]+'"></td></tr>';
      PropTab += '<tr><td><b>Gait</b></td><td>'+properties_species_clicked[2]+'</td><td>'+props_illustrated_species_clicked[2]+'</td></tr></table>';
    html += PropTab
    //show prompt if there is one
    //var prompt = "<p id='prompt' >Given this set of features, to which family does the "+species_clicked+" belong?</p>";
    var prompt = "<p id='prompt' >Given this set of features, to which family does this species belong?</p>";
    html += prompt;
    //display buttons
    var buttons = [];
    if (Array.isArray(trial.button_html)) {
      if (trial.button_html.length == trial.choices.length) {
        buttons = trial.button_html;
      } else {
        console.error('Error in image-button-response plugin. The length of the button_html array does not equal the length of the choices array');
      }
    } else {
      for (var i = 0; i < trial.choices.length; i++) {
        buttons.push(trial.button_html);
      }
    }
    html += '<div id="jspsych-image-button-response-btngroup">';

    for (var i = 0; i < trial.choices.length; i++) {
      var str = buttons[i].replace(/%choice%/g, trial.choices[i]);
      html += '<div class="jspsych-image-button-response-button" style="display: inline-block; margin:'+trial.margin_vertical+' '+trial.margin_horizontal+'" id="jspsych-image-button-response-button-' + i +'" data-choice="'+i+'">'+str+'</div>';
    }
    html += '</div>';

    var btn_placeHolder = "<p id='btnHolder'></p>";
    html += btn_placeHolder;

    display_element.innerHTML = html;

    // start timing
    var start_time = performance.now();

    for (var i = 0; i < trial.choices.length; i++) {
      display_element.querySelector('#jspsych-image-button-response-button-' + i).addEventListener('click', function(e){
        var choice = e.currentTarget.getAttribute('data-choice'); // don't use dataset for jsdom compatibility
        after_response(choice);
      });
    }

    // store response
    var response = {
      rt: null,
      button: null,
    };

    // function to display feedback
    function feedback (feedback_text){document.getElementById('jspsych-image-button-response-btngroup').innerHTML='<span style = font-size:18> '+feedback_text+' </span>'};
    // function to remove prompt
    function removePrompt (){document.getElementById('prompt').innerHTML=''};
   
    // function to handle responses by the subject
    function after_response(choice) {

      // measure rt
      var end_time = performance.now();
      var rt = end_time - start_time;
      response.button = choice;
      response.rt = rt;
      var label_selected = dino_labels_reshuffled_sdl[choice];
      if(label_selected==true_familyName){var choice_correct = "Correct"}else{var choice_correct = "Wrong"}
      response.label_selected = label_selected;
      response.choice_correct = choice_correct;
      // after a valid response, the stimulus will have the CSS class 'responded'
      // which can be used to provide visual feedback that a response was recorded
      display_element.querySelector('#jspsych-image-button-response-stimulus').className += ' responded';

      // disable all the buttons after a response
      var btns = document.querySelectorAll('.jspsych-image-button-response-button button');
      for(var i=0; i<btns.length; i++){
        //btns[i].removeEventListener('click');
        btns[i].setAttribute('disabled', 'disabled');
      }
      
      btn_clickedLabel_index = response.button;
      clickedLabel_name = dino_labels_reshuffled_sdl[btn_clickedLabel_index];
      if(clickedLabel_name==true_familyName){var feedback_text="Your answer is <span style='color:green'><b>correct</b></span>. The family name is <b>"+true_familyName+"</b>."}else{
        feedback_text="Your answer is <span style='color:red'><b>wrong</b></span>. The correct family name is <b>"+true_familyName+"</b>."};
      removePrompt();
      feedback(feedback_text);
      
      function addBtn(){
          document.getElementById('btnHolder').innerHTML = '<input type="button" class="jspsych-btn" id="doneButton2" value="Continue" />';
        };
      addBtn();
      display_element.querySelector('#doneButton2').addEventListener('click', function(e){
        end_trial();
      });
      
      /*if (trial.response_ends_trial) {
        end_trial();
      }*/
    };

    // function to end trial when it is time
    function end_trial() {

      // kill any remaining setTimeout handlers
      jsPsych.pluginAPI.clearAllTimeouts();

      // gather the data to store for the trial
      var trial_data = {
        "rt": response.rt,
        "stimulus": trial.stimulus,
        "button_pressed": response.button,
        "label_selected": response.label_selected,
        "choice_correct": response.choice_correct
      };

      // clear the display
      display_element.innerHTML = '';

      // move on to the next trial
      jsPsych.finishTrial(trial_data);
    };



    // hide image if timing is set
    if (trial.stimulus_duration !== null) {
      jsPsych.pluginAPI.setTimeout(function() {
        display_element.querySelector('#jspsych-image-button-response-stimulus').style.visibility = 'hidden';
      }, trial.stimulus_duration);
    }

    // end trial if time limit is set
    if (trial.trial_duration !== null) {
      jsPsych.pluginAPI.setTimeout(function() {
        end_trial();
      }, trial.trial_duration);
    }

  };

  return plugin;
})();

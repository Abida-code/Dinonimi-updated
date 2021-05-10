/**
 * jspsych-image-button-response
 * Josh de Leeuw
 *
 * plugin for displaying a stimulus and getting a keyboard response
 *
 * documentation: docs.jspsych.org
 *
 **/

jsPsych.plugins["judgmentTest"] = (function() {

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
      choices_conf: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Choices_conf',
        default: undefined,
        array: true,
        description: 'The labels for the buttons.'
      },
      button_html: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Button HTML',
        default: '<button class="jspsych-btn" >%choice%</button>',
        array: true,
        description: 'The html of the button. Can create own style.'
      },
      button_conf_html: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Button HTML',
        default: '<button class="jspsych-btn-confRating" style="color:black">%choice%</button>',
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

    var html = '<div></div>';
    // display number of test item
    html += "<h2>Judgment Knowlege Test (Round "+(count_judgment+1)+"): Question "+((t_judge+1)+(t_judge_block*(n_families*n_subFamilies)))+ " out of " 
      +n_families*n_subFamilies*n_blocks_judgment+"</h2>";

    // display stimulus
    html += '<p><img src="'+trial.stimulus+'" id="jspsych-image-button-response-stimulus" '
    +'style=" height:140px; border:4px solid black; ';
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
    
    //show prompt
    var prompt_judgmentTest = "<p id='prompt_judgmentTest'>Does this dinosaur belong to the family of <b>"
    +names_tobeDisplayed_judgment_curBlock[t_judge]+"</b>?</p>";
    html += prompt_judgmentTest;

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
      html += '<div class="jspsych-image-button-response-button" style="display: inline-block; margin:'
      +trial.margin_vertical+' '
      +trial.margin_horizontal+'" id="jspsych-image-button-response-button-' + i +'" data-choice="'+i+'">'+str+'</div>';
    }
    html += '</div>';

    //display confidence buttons
    html += '<p><div id="prompt_confRating"></div>';
    var buttonsConf = [];
    for (var i = 0; i < trial.choices_conf.length; i++) {
        buttonsConf.push(trial.button_conf_html);
        if(i==0){html += '<p>'};
        html += '<div id="jspsych-image-conf_button-response-btngroup-'
        +i+'" style="display: inline-block;  margin:'+trial.margin_horizontal+'" data-choice-conf="'+i+'"></div>';
    };
    
    display_element.innerHTML = html;

    // start timing
    var start_time = performance.now();

    for (var i = 0; i < trial.choices.length; i++) {
      display_element.querySelector('#jspsych-image-button-response-button-' + i).addEventListener('click', function(e){
        var choice = e.currentTarget.getAttribute('data-choice'); // don't use dataset for jsdom compatibility
        after_response_yn(choice);
      });
    }

    // store response
    var response = {
      rt: null,
      button: null,
      rt_conf: null,
      button_conf: null
    };

    // function to handle responses by the subject
    function after_response_yn(choice) {

      // measure rt
      var end_time = performance.now();
      var rt = end_time - start_time;
      response.button = choice;
      response.rt = rt;
      var choice_label = ['Yes','No'][choice]
      var presentedSpecies = species_sequence_judgment_curBlock[t_judge]
      var presentedLabel = names_tobeDisplayed_judgment_curBlock[t_judge]
      var trueLabel = true_familyNames_judgment_curBlock[t_judge]
      if(trueLabel==presentedLabel){
        if(choice==0){var response_cat="Hit"}else{var response_cat="Miss"}
      }else{
        if(choice==0){var response_cat="FA"}else{var response_cat="CR"}
      }
      response.presentedSpecies = presentedSpecies;
      response.presentedLabel = presentedLabel;
      response.trueLabel = trueLabel;
      response.response_cat = response_cat;
      response.choice_label = choice_label;

      // after a valid response, the stimulus will have the CSS class 'responded'
      // which can be used to provide visual feedback that a response was recorded
      display_element.querySelector('#jspsych-image-button-response-stimulus').className += ' responded';

      // disable all the buttons after a response
      var btns = document.querySelectorAll('.jspsych-image-button-response-button button');
      for(var i=0; i<btns.length; i++){
        //btns[i].removeEventListener('click');
        btns[i].setAttribute('disabled', 'disabled');
        if(i==choice){btns[i].setAttribute('style', 'color: black; border: solid black')};
        
      }

     
      document.getElementById('prompt_confRating').innerHTML='How certain are you in your judgment?';
      
      for (var i = 0; i < trial.choices_conf.length; i++) {
        var str_i = buttonsConf[i].replace(/%choice%/g, trial.choices_conf[i]);
        document.getElementById('jspsych-image-conf_button-response-btngroup-'+i+'').innerHTML=str_i;
      };
      
      // change color of confidence button
      confbtns = document.getElementsByClassName("jspsych-btn-confRating");

      confbtns_colors = ["red","orange","yellow","green"];
      for(i=0;i<confbtns.length;i++){
        confbtns[i].setAttribute('style','background-color:'+confbtns_colors[i]+'');
      };

      for (var i = 0; i < trial.choices_conf.length; i++) {
        display_element.querySelector('#jspsych-image-conf_button-response-btngroup-' + i).addEventListener('click', function(e){
          var choice = e.currentTarget.getAttribute('data-choice-conf'); // don't use dataset for jsdom compatibility
          after_response_confRating(choice);
        });
      }
      
    };

    function after_response_confRating(choice) {
      
      var end_time_conf = performance.now();
      var rt_conf = end_time_conf - start_time;
      response.button_conf = choice;
      var chosen_confRating = trial.choices_conf[choice]
      response.chosen_confRating = chosen_confRating;
      response.rt_conf = rt_conf;

      if (trial.response_ends_trial) {
        end_trial();
      }
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
        "choice_label": response.choice_label,
        "presentedSpecies": response.presentedSpecies,
        "presentedLabel": response.presentedLabel,
        "trueLabel": response.trueLabel,
        "response_cat": response.response_cat,
        "button_pressed_conf": response.button_conf,
        "chosen_confRating": response.chosen_confRating,
        "rt_conf": response.rt_conf
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

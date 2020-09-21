// var testArray = [
// {index: 0, string: 'string75', value: 95, weight: 75},
// {index: 1, string: 'string25', value: 5, weight: 25}
// ]

// var Index = []
// var Weight = []

// for (i=0;i<testArray.length;i++) {
//     Index[i] = testArray[i].index;
//     Weight[i] = testArray[i].weight;
// }

// class WeightedSampler {
//   constructor(elements, weights) {
//     this.total = 0;
//     this.elements = Array.from(elements);
//     this.cweights = weights.map(weight => this.total += weight);
//   }
//   get() {
//     let random = Math.random() * this.total;
//     return this.elements.find((element, index) => random < this.cweights[index]);
//   }
// }

// const sampler8 = new WeightedSampler(Index, Weight);
// var randomArray = Array.apply(null, Array(20)).map(() => sampler8.get());

// var finalArray = []

// for  (i=0;i<randomArray.length;i++)  {
// 	finalArray[i] = testArray[randomArray[i]];
// }

// finalArray[0] = testArray[randomArray[0]]

// finalArray[1] = testArray[randomArray[1]]


// AllChamps = _.sample(ChampArray, 10)

// for(i = 0; i <= num_comps; i++){
// AllChamps = _.sample(ChampArray, 10)
// BlueTeam[i] = AllChamps.slice(0,5);
// RedTeam[i] = AllChamps.slice(5,10);
// }

// exp.PresentItems = []

// exp.PresentItems = LoLStimuli

// exp.PresentItems[0] = {Blue: BlueTeam[0], Red: RedTeam[0]}

// for(i = 0; i < num_comps; i++){
// 	exp.PresentItems[i] = {Blue: BlueTeam[i], Red: RedTeam[i]}
// 	}

function make_slides(f) {
  var   slides = {};

  slides.i0 = slide({
     name : "i0",
     start: function() {
      exp.startT = Date.now();
     }
  });

  slides.instructions = slide({
    name : "instructions",
    button : function() {
      exp.go(); //use exp.go() if and only if there is no "present" data.
    }
  });

  // slides.single_trial = slide({
  //   name: "single_trial",
  //   start: function() {
  //     $(".err").hide();
  //     $(".display_condition").html("You are in " + exp.condition + ".");
  //   },
  //   button : function() {
  //     response = $("#text_response").val();
  //     if (response.length == 0) {
  //       $(".err").show();
  //     } else {
  //       exp.data_trials.push({
  //         "trial_type" : "single_trial",
  //         "response" : response
  //       });
  //       exp.go(); //make sure this is at the *end*, after you log your data
  //     }
  //   },
  // });

  slides.main_questions = slide({
    name : "main_questions",

    /* trial information for this block
     (the variable 'stim' will change between each of these values,
      and for each of these, present_handle will be run.) */
    present : LoLStimuli,
        // "bins" : [
        //   {
        //     "min" : 0,
        //     "max" : 10
        //   },
        //   {
        //     "min" : 10,
        //     "max" : 20
        //   },
        //   {
        //     "min" : 20,
        //     "max" : 30
        //   },
        //   {
        //     "min" : 30,
        //     "max" : 40
        //   },
        //   {
        //     "min" : 40,
        //     "max" : 50
        //   },
        //   {
        //     "min" : 50,
        //     "max" : 60
        //   }
        // ],
        // "question": "How tall is tall?",
    // present : [exp.PresentItems[0]],
    // present : _.sample(exp.PresentItems),

    //this gets run only at the beginning of the block
    present_handle : function(stim) {
      $(".err_radio1").hide();
      $(".err_radio2").hide();
      $(".err_slider").hide();
      // $(".showButton").show();

      this.stim = stim; //I like to store this information in the slide so I can record it later.

      $(".BlueTeamLabel").html("Blue Team");
      $(".RedTeamLabel").html("Red Team");

      $(".Blue1Name").html(stim.Blue[0].Name)
      $(".Blue2Name").html(stim.Blue[1].Name)
      $(".Blue3Name").html(stim.Blue[2].Name)
      $(".Blue4Name").html(stim.Blue[3].Name)
      $(".Blue5Name").html(stim.Blue[4].Name)


      $(".Red1Name").html(stim.Red[0].Name)
      $(".Red2Name").html(stim.Red[1].Name)
      $(".Red3Name").html(stim.Red[2].Name)
      $(".Red4Name").html(stim.Red[3].Name)
      $(".Red5Name").html(stim.Red[4].Name)

      $("#Blue1Image").html("<img src =\"" + stim.Blue[0].Image + "\" alt=\"Blue1\" id=\"Blue1Image\"></img>")
      $("#Blue2Image").html("<img src =\"" + stim.Blue[1].Image + "\" alt=\"Blue2\" id=\"Blue2Image\"></img>")
      $("#Blue3Image").html("<img src =\"" + stim.Blue[2].Image + "\" alt=\"Blue3\" id=\"Blue3Image\"></img>")
      $("#Blue4Image").html("<img src =\"" + stim.Blue[3].Image + "\" alt=\"Blue4\" id=\"Blue4Image\"></img>")
      $("#Blue5Image").html("<img src =\"" + stim.Blue[4].Image + "\" alt=\"Blue5\" id=\"Blue5Image\"></img>")

      // $("#Red1Image").html("<img src =\"" + stim.Red[0].Image + "\" alt=\"Red1\" id=\"Red1Image\"></img>")
      // $("#Red2Image").html("<img src =\"" + stim.Red[1].Image + "\" alt=\"Red2\" id=\"Red1Image\"></img>")
      // $("#Red3Image").html("<img src =\"" + stim.Red[2].Image + "\" alt=\"Red3\" id=\"Red1Image\"></img>")
      // $("#Red4Image").html("<img src =\"" + stim.Red[3].Image + "\" alt=\"Red4\" id=\"Red1Image\"></img>")
      // $("#Red5Image").html("<img src =\"" + stim.Red[4].Image + "\" alt=\"Red5\" id=\"Red1Image\"></img>")

      $(".EarlyStatement").html(stim.EarlyStatement)
      $(".LateStatement").html(stim.LateStatement)

	this.sentence_types = ["Early", "Late", "General"];
      var sentences = {
        "Early": "...excel in the laning phase?",
        "Late": "...excel in the late game?",
        "General": "...excel in general?",
      };

  this.sentence_types2 = ["Early"];
      var sentences2 = {
        "Early": "...excelled in the laning phase?",
      };

  this.sentence_types3 = [ "Late"];
      var sentences3 = {
        "Late": "...excelled in the late game?",
      };

      this.n_sliders_Blue = this.sentence_types.length;
      $(".slider_row_Blue").remove();
      for (var i=0; i<this.n_sliders_Blue; i++) {
        var sentence_type = this.sentence_types[i];
        var sentence = sentences[sentence_type];
        $("#multi_slider_table_Blue").append('<tr class="slider_row_Blue"><td class="slider_target" id="objectBlue' + i + '">' + sentence + '</td><td colspan="2"><div id="sliderBlue' + i + '" class="slider">-------[ ]--------</div></td></tr>');
        utils.match_row_height("#multi_slider_table_Blue", ".slider_target");
      }

      this.n_sliders_Speaker_Early = this.sentence_types2.length;
      $(".slider_row_Speaker_Early").remove();
      for (var i=0; i<this.n_sliders_Speaker_Early; i++) {
        var sentence_type2 = this.sentence_types2[i];
        var sentence = sentences2[sentence_type2];
        $("#multi_slider_table_Speaker_Early").append('<tr class="slider_row_Speaker_Early"><td class="slider_target" id="objectSpeaker_Early' + i + '">' + sentence + '</td><td colspan="2"><div id="sliderSpeaker_Early' + i + '" class="slider">-------[ ]--------</div></td></tr>');
        utils.match_row_height("#multi_slider_table_Speaker_Early", ".slider_target");
      }

      this.n_sliders_Speaker_Late = this.sentence_types3.length;
      $(".slider_row_Speaker_Late").remove();
      for (var i=0; i<this.n_sliders_Speaker_Late; i++) {
        var sentence_type3 = this.sentence_types3[i];
        var sentence = sentences3[sentence_type3];
        $("#multi_slider_table_Speaker_Late").append('<tr class="slider_row_Speaker_Late"><td class="slider_target" id="objectSpeaker_Late' + i + '">' + sentence + '</td><td colspan="2"><div id="sliderSpeaker_Late' + i + '" class="slider">-------[ ]--------</div></td></tr>');
        utils.match_row_height("#multi_slider_table_Speaker_Late", ".slider_target");
      }


      // this.init_sliders(this.preferences);
      this.init_sliders_Blue(this.sentence_types);
      this.init_sliders_Speaker_Early(this.sentence_types2);
      this.init_sliders_Speaker_Late(this.sentence_types3)
      exp.sliderPostBlue = [];
      exp.sliderPostSpeaker_Early = [];
      exp.sliderPostSpeaker_Late = [];

      $('input[type=radio]').attr('checked', false); //for radio button response
      // hide stuff
      $(".err_radio1").hide();
      $(".err_radio2").hide();
      $(".err_slider").hide();
      $(".err_slider2").hide();
      $(".err_slider3").hide();
      $(".hidden1").hide();
      $(".hidden2").hide();
      $(".hidden3").hide();
      $(".hidden4").hide();

      this.trial_num++;
    },

    showButton : function() {
      if ($('input[name=ChampFam]:checked').size() <= 0) {
        $(".err_radio1").show();
      } else {
        this.log_responses();
        $(".showButton").hide();
        // $(".radio").show();
        $(".hidden1").show();
        $(".err_radio1").hide();                   
        }

        /* use _stream.apply(this); if and only if there is
        "present" data. (and only *after* responses are logged) */
        // _stream.apply(this);
      //response = $("#testFreeResponse").val();
    },

	init_sliders_Blue : function() {
      for (var i=0; i<this.sentence_types.length; i++) {
         utils.make_slider("#sliderBlue" + i, this.make_slider_callback_Blue(i));
      }
    },
    init_sliders_Speaker_Early : function() {
      for (var i=0; i<this.sentence_types2.length; i++) {
         utils.make_slider("#sliderSpeaker_Early" + i, this.make_slider_callback_Speaker_Early(i));
      }
    },
  init_sliders_Speaker_Late : function() {
      for (var i=0; i<this.sentence_types3.length; i++) {
         utils.make_slider("#sliderSpeaker_Late" + i, this.make_slider_callback_Speaker_Late(i));
      }
    },
    make_slider_callback_Blue : function(i) {
      return function(event, ui) {
        exp.sliderPostBlue[i] = ui.value;
      };
    },
    make_slider_callback_Speaker_Early : function(i) {
      return function(event, ui) {
        exp.sliderPostSpeaker_Early[i] = ui.value;
      };
    },
    make_slider_callback_Speaker_Late : function(i) {
      return function(event, ui) {
        exp.sliderPostSpeaker_Late[i] = ui.value;
      };
    },

    showButton2: function(){

      if ($('input[name=TVJearly]:checked').size() <= 0 || $('input[name=TVJlate]:checked').size() <= 0 ) {
        $(".err_radio2").show();
      } 
      else{  
        this.log_responses();
        $(".showButton2").hide();
        $(".hidden2").show();
        $(".err_radio2").hide();
      }
    },

    showButton3 : function(){
      var ok_to_go_on = true
      for (var i=0; i<this.n_sliders_Blue; i++) {
        if (exp.sliderPostBlue[i]==undefined){
          ok_to_go_on = false
        }
      }
      if(ok_to_go_on){
        this.log_responses();
        $(".showButton3").hide();
        $(".hidden3").show();
        $(".err_slider").hide();
      }
      else{
        $(".err_slider").show();
      }
    },

    showButton4 : function(){
      var ok_to_go_on = true
      for (var i=0; i<this.n_sliders_Speaker_Early; i++) {
        if (exp.sliderPostSpeaker_Early[i]==undefined){
          ok_to_go_on = false
        }
      }
      if(ok_to_go_on){
        this.log_responses();
        $(".showButton4").hide();
        $(".hidden4").show();
        $(".err_slider2").hide();
      }
      else{
        $(".err_slider2").show();
      }
    },


    button : function(){
      	var ok_to_go_on = true
      	for (var i=0; i<this.n_sliders_Speaker_Late; i++) {
        if (exp.sliderPostSpeaker_Late[i]==undefined){
          ok_to_go_on = false
        }
      }
      if(ok_to_go_on){
      	var end_time = Date.now();
        this.time_spent = end_time - this.start_time;
        this.log_responses();
        $(".showButton").show();
        _stream.apply(this); //make sure this is at the *end*, after you log your data
      }
      	else {
           $(".err_slider3").show();         
      }
    },

    log_responses : function() {
      exp.data_trials.push({
        "blueTopIndex" : this.stim.Blue[0].Index,
        "blueTopName" : this.stim.Blue[0].Name,
        "blueTopPlay" : this.stim.Blue[0].Top,
        "blueJungleIndex" : this.stim.Blue[1].Index,
        "blueJungleName" : this.stim.Blue[1].Name,
        "blueJunglePlay" : this.stim.Blue[1].Jungle,
        "blueMidIndex" : this.stim.Blue[2].Index,
        "blueMidName" : this.stim.Blue[2].Name,
        "blueMidPlay" : this.stim.Blue[2].Mid,
        "blueBotIndex" : this.stim.Blue[3].Index,
        "blueBotName" : this.stim.Blue[3].Name,
        "blueBotPlay" : this.stim.Blue[3].Bot,
        "blueSupportIndex" : this.stim.Blue[4].Index,
        "blueSupportName" : this.stim.Blue[4].Name,
        "blueSupportPlay" : this.stim.Blue[4].Support,
        
        "blue_excel_early" : exp.sliderPostBlue[0],
        "blue_excel_late" : exp.sliderPostBlue[1],
        "blue_excel_general" : exp.sliderPostBlue[2],
        // "realistic_blue" : $('input[name=realistic_blue]:checked').val(), //if using radio buttons
        "speaker_early_index" : this.stim.EarlyStatementIndex,
        "speaker_late_index" : this.stim.LateStatementIndex,
        "speaker_excel_early": exp.sliderPostSpeaker_Early[0],
        "speaker_excel_late": exp.sliderPostSpeaker_Late[0],
        "ChampFam" : $('input[name=ChampFam]:checked').val(), //if using radio buttons
        "TVJearly" : $('input[name=TVJearly]:checked').val(), //if using radio buttons
        "TVJlate" : $('input[name=TVJlate]:checked').val() //if using radio buttons
      });
    }
  });

slides.expertise = slide({
  name: "expertise",

  present: [{"emptiness": "dear god please work"}],

present_handle : function(stim){
$(".err_mega").hide();

this.sentence_types = ["Early", "Late"];
      var sentences = {
        "Early": "...excel in the laning phase?",
        "Late": "...excel in the late game?",
      };

      this.n_sliders_Expertise = this.sentence_types.length;
      $(".slider_row_Expertise").remove();
      for (var i=0; i<this.n_sliders_Expertise; i++) {
        var sentence_type = this.sentence_types[i];
        var sentence = sentences[sentence_type];
        $("#multi_slider_table_Expertise").append('<tr class="slider_row_Expertise"><td class="slider_target" id="objectExpertise' + i + '">' + sentence + '</td><td colspan="2"><div id="sliderExpertise' + i + '" class="slider">-------[ ]--------</div></td></tr>');
        utils.match_row_height("#multi_slider_table_Expertise", ".slider_target");
      }

      this.init_sliders_Expertise(this.sentence_types);
      exp.sliderPostExpertise = [];
  },

  init_sliders_Expertise : function() {
      for (var i=0; i<this.sentence_types.length; i++) {
         utils.make_slider("#sliderExpertise" + i, this.make_slider_callback_Expertise(i));
      }
    },

  make_slider_callback_Expertise : function(i) {
      return function(event, ui) {
        exp.sliderPostExpertise[i] = ui.value;
      };
    },

    megaButton : function(){
        var ok_to_go_on = true
        for (var i=0; i<this.n_sliders_Expertise; i++) {
        if (exp.sliderPostExpertise[i]==undefined){
          ok_to_go_on = false
        }
      }
      if(ok_to_go_on && $("#year_started").val()!="-1" && $("#league_level").val()!= "-1" && $("#hours_total").val()!= "-1"
        && $("#games_weekly").val()!=undefined && $("#champfam_total").val()!="-1"){
        var end_time = Date.now();
        // this.time_spent = end_time - this.start_time;
        // this.log_responses();
        _stream.apply(this); //make sure this is at the *end*, after you log your data
      }
        else {
           $(".err_mega").show();         
      }
    },

  log_responses : function() {
    exp.subj_expertise = {
    "year_started" : $("#year_started").val(),
    "league_level" : $("#league_level").val(),
    "hours_total" : $("#hours_total").val(),
    "games_weekly" : $("#games_weekly").val(),
    "champfam_total" : $("champfam_total").val(),
    "hierarch_early" : exp.sliderPostExpertise[0],
    "hierarch_late" : exp.sliderPostExpertise[1]
  }
    exp.go(); //use exp.go() if and only if there is no "present" data.
  } 

})

// slides.prior_knowledge = slide({
//     name : "prior_knowledge",

//      trial information for this block
//      (the variable 'stim' will change between each of these values,
//       and for each of these, present_handle will be run.) 
//     present : exp.PresentItems2,

    // trial_num: 0,

    //this gets run only at the beginning of the block
    // present_handle : function(stim) {

    //   $('input[type=radio]').attr('checked', false); //for radio button response
    //   // hide stuff
    //   $(".err").hide();

    //   this.stim = stim; //I like to store this information in the slide so I can record it later.


    //   // $(".CapFact").html(stim.Fact.charAt(0).toUpperCase() + stim.Fact.substring(1)+".");
    //   $(".QuestFact").html(stim.Fact + "?")
    //   $(".Animal").html(stim.Animal)
    //   "<img src=" + stim.Image + "alt=\"Animal\" id=\"AnimalPic\"></img>"
    //   $("#AnimalExemplar2").html("<img src =\"" + stim.Image + "\" alt=\"Animal\" id=\"AnimalPic\"></img>")
    //         this.start_time = Date.now()
    //   // this.init_sliders();
    //   // exp.sliderPost = null; //erase current slider value
    //   this.trial_num++;
    // },

    //     button : function() {
    //   var end_time = Date.now();
    //   //response = $("#testFreeResponse").val();
    //   if ($('input[type=radio]:checked').size() == 0) {
    //     $(".err").show();
    //   } else {
    //     this.time_spent = end_time - this.start_time;
    //     this.log_responses();
    //     _stream.apply(this); //make sure this is at the *end*, after you log your data
    //   }
    // },

  //   log_responses : function() {
  //     exp.data_trials.push({
  //       "trial_type" : "prior_knowledge",
  //       "response" : $('input[type=radio]:checked').val(), //if using radio buttons
  //       // "Animal" : stim.Animal,
  //       "Format" : this.stim.Format,
  //       "TraitType" : this.stim.TraitType,
  //       "Fact" : this.stim.Fact
  //     });
  //   }
  // });

  slides.subj_info =  slide({
    name : "subj_info",
    
    submit : function(e){
      //if (e.preventDefault) e.preventDefault(); // I don't know what this means.
      exp.subj_data = {
        language : $("#language").val(),
        // enjoyment : $("#enjoyment").val(),
        // asses : $('input[name="assess"]:checked').val(),
        age : $("#age").val(),
        gender : $("#gender").val(),
        education : $("#education").val(),
        comments : $("#comments").val()
        // problems: $("#problems").val(),
        // fairprice: $("#fairprice").val()
      };
      exp.go(); //use exp.go() if and only if there is no "present" data.
    }
  });

  slides.thanks = slide({
    name : "thanks",
    start : function() {
      exp.data= {
          "trials" : exp.data_trials,
          "catch_trials" : exp.catch_trials,
          "system" : exp.system,
          "condition" : exp.condition,
          "subject_information" : exp.subj_data,
          "subject_expertise" : exp.subj_expertise,
          "time_in_minutes" : (Date.now() - exp.startT)/60000
      };
      setTimeout(function() {turk.submit(exp.data);}, 1000);
    }
  });

  return slides;
}

/// init ///
function init() {
  exp.trials = [];
  exp.catch_trials = [];
  exp.condition = _.sample(["CONDITION 1", "condition 2"]); //can randomize between subject conditions here
  exp.system = {
      Browser : BrowserDetect.browser,
      OS : BrowserDetect.OS,
      screenH: screen.height,
      screenUH: exp.height,
      screenW: screen.width,
      screenUW: exp.width
    };
  //blocks of the experiment:
  exp.structure=["i0", "instructions", "main_questions", "expertise", 'subj_info', 'thanks'];

  exp.data_trials = [];
  //make corresponding slides:
  exp.slides = make_slides(exp);

  exp.nQs = utils.get_exp_length(); //this does not work if there are stacks of stims (but does work for an experiment with this structure)
                    //relies on structure and slides being defined

  $('.slide').hide(); //hide everything

  //make sure turkers have accepted HIT (or you're not in mturk)
  $("#start_button").click(function() {
    if (turk.previewMode) {
      $("#mustaccept").show();
    } else {
      $("#start_button").click(function() {$("#mustaccept").show();});
      exp.go();
    }
  });

  exp.go(); //show first slide
}
